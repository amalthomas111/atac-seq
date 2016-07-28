#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  Author: Amal Thomas

import os, sys, json, time, subprocess

if len(sys.argv)!=3:
	print("Usage: python3 preprocessing.py <jsonfile> <inputfilename>\nExiting!!!!\n")
	exit(0)
if not os.path.exists(sys.argv[1]):
	print("Error: Input jsonfile not found \nExiting!!\n")
	exit(0)

INPUT = sys.argv[2].strip()
def check_directory(param):
	#check input file and folders
	if not os.path.isdir(param["folders"]["rawreads"]) :
		print("rawreads directory not found\nExiting!!!!\n")
		exit(0)
	if not os.path.exists(param["folders"]["rawreads"]+"/"+sys.argv[2]+"_1.fastq") or \
	not os.path.exists(param["folders"]["rawreads"]+"/"+sys.argv[2]+"_2.fastq") :
		print("Input file not found\nExiting!! \n")
		exit(0)
	if not os.path.exists(param["folders"]["reads"]):
		os.makedirs(param["folders"]["reads"])
	if not os.path.exists(param["folders"]["bams"]):
		os.makedirs(param["folders"]["bams"])
	if not os.path.exists(param["folders"]["beds"]):
		os.makedirs(param["folders"]["beds"])
	if not os.path.exists(param["folders"]["tracks"]):
		os.makedirs(param["folders"]["tracks"])
	if not os.path.exists(param["folders"]["peaks"]):
		os.makedirs(param["folders"]["peaks"])
	if not os.path.exists(param["folders"]["quality"]):
		os.makedirs(param["folders"]["quality"])

def quality_checks(param):
	#trim adaptors
	trim_adaptors = param["programs"]["trim_galore"]+" --paired --nextera "+ \
	param["folders"]["rawreads"]+"/"+ INPUT +"_1.fastq "+ \
	param["folders"]["rawreads"]+"/"+ INPUT +"_2.fastq "+ \
	"--output_dir "+ param["folders"]["reads"]
	
	#fastqc quality analyzis
	fastqc = param["programs"]["fastqc"]+ " " +\
	param["folders"]["reads"]+ "/" + INPUT + "_1_val_1.fq " +\
	param["folders"]["reads"]+ "/" + INPUT + "_2_val_2.fq -o "+ \
	param["folders"]["quality"]
	
	#renaming files
	rename1 = "mv "+param["folders"]["reads"]+ "/" + INPUT + "_1_val_1.fq "+\
	param["folders"]["reads"]+ "/" + INPUT +"_1_trimmed.fastq"
	rename2 = "mv "+param["folders"]["reads"]+ "/" + INPUT + "_2_val_2.fq "+\
	param["folders"]["reads"]+ "/" + INPUT +"_2_trimmed.fastq"
	
	print(trim_adaptors,fastqc,rename1,rename2,sep="\n")

	os.system(trim_adaptors)
	os.system(fastqc)
	os.system(rename1)
	os.system(rename2)

def mapping(param):
	#bowtie alignment
	bowtie = param["programs"]["bowtie2"]+" -p 8 -X 2000 --fr --no-discordant --no-mixed --minins 38 "+\
	"--met-file "+ param["folders"]["quality"]+ "/"+ INPUT + ".alignmetrics.txt "+\
	"-x  " + param["programs"]["bowtie_index"] + \
	" -1 " + param["folders"]["reads"]+ "/" + INPUT + "_1_trimmed.fastq"+ \
	" -2 " + param["folders"]["reads"]+ "/" + INPUT + "_2_trimmed.fastq |"+ \
	param["programs"]["samtools"]+" view -bS - | "+param["programs"]["samtools"]+ " sort - "+\
	param["folders"]["bams"]+"/"+ INPUT + ".sorted"
	print(bowtie)
	os.system(bowtie)
def quality_controls(param):
	#filter unwanted reads
	remove_junk_chrom = "samtools view -h "+\
	param["folders"]["bams"]+"/"+ INPUT + ".sorted.bam |"+\
	'''awk 'substr($0,1,1) == "@"|| (length($3)<=5 && $3!="chrM" && $3 != "*")'''+\
	'''{{print $0}}'| '''+ param["programs"]["samtools"]+ ''' view -bS - > ''' +\
	param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.bam "
	
	#remove duplicate reads
	remove_duplicate = "java -Xmx2g -jar "+param["programs"]["markduplicates"]+" "+\
	"INPUT="+param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.bam "+\
	"METRICS_FILE="+param["folders"]["bams"]+"/" + INPUT + ".markdup.metrics " +\
	"OUTPUT="+param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.nodup.bam " +\
	"REMOVE_DUPLICATES=true "+\
	"ASSUME_SORTED=true "
	
	#index bam file
	index_bam = param["programs"]["samtools"]+" index "+\
	param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.nodup.bam "
	
	#mark duplicate reads
	mark_duplicate = "java -Xmx2g -jar "+param["programs"]["markduplicates"]+" "+\
	"INPUT="+param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.bam "+\
	"METRICS_FILE="+param["folders"]["bams"]+"/" + INPUT + ".markdup.metrics " +\
	"OUTPUT="+param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.markdup.bam " +\
	"REMOVE_DUPLICATES=false "+\
	"ASSUME_SORTED=true "
	
	#remove low quality reads
	remove_lowquality = param["programs"]["samtools"]+"  view -b -h -q 30 "+\
	param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.nodup.bam  > " +\
	param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.nodup.filt.bam"
	
	#qualimap stats
	qualimap = param["programs"]["qualimap"]+ " bamqc -bam "+\
	param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.nodup.bam "+\
	"--outdir "+param["folders"]["quality"]+" --outfile "+\
	INPUT  + ".filt.sorted.chr.nodup.qualimap.pdf"
	
	print(remove_junk_chrom,index_bam,remove_duplicate,mark_duplicate,remove_lowquality,qualimap,sep="\n")
	
	os.system(remove_junk_chrom)
	os.system(index_bam)
	os.system(remove_duplicate)
	os.system(mark_duplicate)
	os.system(remove_lowquality)
	os.system(qualimap)
def peak_calling(param):
	#shift cordinates
	temp_bed = param["programs"]["bamtobed"]+ " -i "+ \
	param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.nodup.filt.bam > "+ INPUT + ".temp.bed"
	insert_temp = param["programs"]["samtools"]+" view "+\
	param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.nodup.filt.bam |"+\
	''' awk 'BEGIN{OFS="\t"}{print $1,$9}' >  '''+\
	INPUT +".insert.temp"
	shift = "paste "+ INPUT + ".temp.bed " + INPUT +".insert.temp |"+\
	'''awk 'BEGIN{OFS="\t"}{split($4,name,"/");if(name[1]== $7 && $6 == "+" ) '''+\
	'''{print $1,$2+4,$3,$4,$8,$6}else if(name[1]== $7 && $6 == "-" ){print $1,$2,$3-5,$4,$8,$6} '''+\
	''' else {print "ERROR",name[1],$7}}' > ''' + \
	param["folders"]["beds"]+"/"+ INPUT + ".sorted.chr.nodup.filt.shift.bed"
	remove_temp = "rm " + INPUT + ".temp.bed "+ INPUT +".insert.temp"
	
	#macs peak calling
	macs = "export "+ param["env"]["macs_python_env"]+"; "+\
	param["programs"]["macs2"]+ " callpeak --gsize hs -f BED "+\
	"--keep-dup all  --call-summits  --shift -100 --extsize 200 "+\
	"--nomodel --nolambda --verbose 3 "+\
	"--treatment "+ param["folders"]["beds"]+ "/"+ INPUT + ".sorted.chr.nodup.filt.shift.bed "+  \
	"--outdir "+ param["folders"]["peaks"]+ \
	" --name "+ INPUT 

	#remove blacklisted regions from peaks
	remove_blacklist= "intersectBed -v -a "+ param["folders"]["peaks"]+"/"+ INPUT + \
	"_peaks.narrowPeak -b "+param["files"]["blacklist_regions"]+ " >  "+\
	param["folders"]["peaks"]+"/"+ INPUT +"_peaks.narrowPeak.blacklistcleared"
	
	print(temp_bed, insert_temp,shift , remove_temp,
	 macs,remove_blacklist,sep="\n")
	
	os.system(temp_bed)
	os.system(insert_temp)
	os.system(shift)
	os.system(remove_temp)
	#os.system(env)
	os.system(macs)
	os.system(remove_blacklist)
def visualization(param):
	#create bigbed
	temp_track = ''' awk '{{OFS="\t"; print $1, $2, $3}}' '''+\
	param["folders"]["peaks"]+"/"+ INPUT +"_peaks.narrowPeak.blacklistcleared | "+\
	 "sort -k 1,1 -k 2,2n > "+ param["folders"]["tracks"]+"/"+ INPUT + ".temp"
	create_bb = param["programs"]["bedToBigBed"] + " "+param["folders"]["tracks"]+"/"+ INPUT + ".temp "+ \
	param["files"]["hg19.sizes.txt"]+ "  "+param["folders"]["tracks"]+"/"+ INPUT + ".bb"
	
	#create center bed
	center = "cat "+ param["folders"]["beds"]+ "/"+ INPUT + ".sorted.chr.nodup.filt.shift.bed |"+\
	'''awk 'BEGIN{OFS="\t"}{if($6 == "-") $2=$3-1; print $1, $2, $2+1, $4, $5, $6}'|'''+\
	'''sort -k 1,1 -k2,2n > '''+\
	param["folders"]["beds"]+ "/"+ INPUT + ".sorted.chr.nodup.filt.shift.center.bed "
	#extend reads around cut site center
	extend = "bedtools slop -i "+ param["folders"]["beds"]+ "/"+ INPUT +\
	 ".sorted.chr.nodup.filt.shift.center.bed -g "+ param["files"]["hg19.sizes.txt"]+ " -b 15 > "\
	+ param["folders"]["beds"]+ "/"+ INPUT + ".15bpshifted.bed" 
	#libsize
	libsize = "librarySize=$("+ param["programs"]["samtools"]+" view -c -F 4 "+\
	param["folders"]["bams"]+"/"+ INPUT + ".sorted.chr.nodup.filt.bam); "+''' expr="1000000 / $librarySize"; '''+\
	"scaling_factor=$(echo $expr | bc -l); touch "+ param["folders"]["quality"]+"/"+ INPUT +\
	"_lib.size.$librarySize;"+ param["programs"]["bedtools"] + " genomecov -i "+ param["folders"]["beds"]+ "/"+ INPUT +\
	 ".15bpshifted.bed -g "+ param["files"]["hg19.sizes.txt"]+ " -bg -scale $scalingfactor > "+\
	 param["folders"]["beds"]+ "/"+ INPUT +".bedgraph"
	 
	#make bigwigh
	bw = param["programs"]["bedGraphToBigWig"]+ " "+param["folders"]["beds"]+ "/"+ INPUT +".bedgraph "\
	+ param["files"]["hg19.sizes.txt"]+"  "+param["folders"]["tracks"]+"/"+ INPUT + ".bw"
	
	remove_temp = "rm "+param["folders"]["tracks"]+"/"+ INPUT + ".temp"
	print(temp_track,create_bb,center,extend, libsize, bw, remove_temp, sep="\n")
	os.system(temp_track)
	os.system(create_bb)
	os.system(center)
	os.system(extend)
	os.system(libsize)
	#os.system(bedgraph)
	#subprocess.Popen(libsize+bedgraph)
	os.system(bw)
	os.system(remove_temp)
def main():
	with open("preprocessing_setup.json") as f:
		param = json.load(f)
	check_directory(param)
	print("\n",time.strftime("%d-%m-%Y %H:%M:%S", time.localtime()),"\n","### START ###",sep='')
	quality_checks(param)
	mapping(param)
	quality_controls(param)
	peak_calling(param)
	visualization(param)
	print("\n### DONE ###","\n",time.strftime("%d-%m-%Y %H:%M:%S", time.localtime()),sep='')
if __name__=="__main__":
    main()
