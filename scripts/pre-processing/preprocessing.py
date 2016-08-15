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
	print("Error: Input jsonfile not found \n"+
	"Usage:python3 preprocessing.py <jsonfile> <inputfilename> \nExiting!!\n")
	exit(0)

INPUT = sys.argv[2].strip()

with open("preprocessing_setup.json") as f:
	param = json.load(f)
#directories
rawreadsDir = os.path.abspath(param["folders"]["rawreads"])
readsDir = os.path.abspath(param["folders"]["reads"])
qualityDir = os.path.abspath(param["folders"]["quality"])
bamsDir = os.path.abspath(param["folders"]["bams"])
bedsDir = os.path.abspath(param["folders"]["beds"])
peaksDir = os.path.abspath(param["folders"]["peaks"])
tracksDir = os.path.abspath(param["folders"]["tracks"])

print(readsDir)
def check_directory():
	#check input file and folders
	if not os.path.isdir(rawreadsDir) :
		print("rawreads directory not found\nExiting!!!!\n")
		exit(0)
	fastq1 = os.path.join(rawreadsDir, INPUT + "_R1.fastq")
	fastq2 = os.path.join(rawreadsDir, INPUT + "_R2.fastq")
	if not os.path.exists(fastq1) or not os.path.exists(fastq2) :
		print("Input fastq file(s) not found\nExiting!! \n")
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


def quality_checks():
	#trim adaptors
	trim_parameters = {
	'trim_galore' :  param["programs"]["trim_galore"],
	'fastq1' : os.path.join(rawreadsDir, INPUT + "_R1.fastq"),
	'fastq2' : os.path.join(rawreadsDir, INPUT + "_R2.fastq"),
	'trim_parameter' : "--paired --nextera",
	'destinationDir' : readsDir
	}
	trimCMD = "{trim_galore} {trim_parameter} {fastq1} {fastq2} --output_dir {destinationDir}".format(**trim_parameters)
	print(trimCMD)
	os.system(trimCMD)

	#fastqc quality analyzis
	fastqc_quality = {
	'fastqc' :  param["programs"]["fastqc"],
	'trimmed_file1' : os.path.join(readsDir, INPUT + "_R1_val_1.fq"),
	'trimmed_file2' : os.path.join(readsDir, INPUT + "_R2_val_2.fq"),
	'destinationDir' :qualityDir
	}
	fastqcCMD = "{fastqc} {trimmed_file1} {trimmed_file2} -o {destinationDir}".format(**fastqc_quality)
	print(fastqcCMD)
	os.system(fastqcCMD)
	
	#renaming files
	file1 = os.path.join(readsDir, INPUT + "_R1_val_1.fq")
	file2 = os.path.join(readsDir, INPUT + "_R2_val_2.fq")
	file1_new = os.path.join(readsDir, INPUT + "_1_trimmed.fastq")
	file2_new = os.path.join(readsDir, INPUT + "_2_trimmed.fastq")

	os.rename(file1, file1_new)
	os.rename(file2, file2_new)

def mapping():
	#bowtie alignment
	bowtie = {
	'bowtie2' :  param["programs"]["bowtie2"],
	'mapping_parameters': "-p 8 -X 2000 --fr --no-discordant --no-mixed --minins 38",
	'bowtie_index': param["programs"]["bowtie_index"],
	'mate1' : os.path.join(readsDir, INPUT + "_1_trimmed.fastq"),
	'mate2' : os.path.join(readsDir, INPUT + "_2_trimmed.fastq"),
	'metric_file' : os.path.join(qualityDir, INPUT + ".alignmetrics.txt"),
	'samtools' : param["programs"]["samtools"],
	'outputbam' : os.path.join(bamsDir, INPUT + ".sorted")
	}
	bowtieCMD = "{bowtie2} {mapping_parameters} --met-file {metric_file} -x {bowtie_index} \
	-1 {mate1} -2 {mate2} | {samtools} view -bS - | {samtools} sort - {outputbam}".format(**bowtie)
	 
	print(bowtieCMD)
	os.system(bowtieCMD)

def remove_junk_chromosome():
	#function to remove unwanted chromosomes
	remove_chroms = {
	'sortedbam' :  os.path.join(bamsDir, INPUT + ".sorted.bam"),
	'outputbam': os.path.join(bamsDir, INPUT + ".sorted.chr.bam"),
	'tempsam' : "tempsam_"+ INPUT,
	'samtools' : param["programs"]["samtools"]
	}
	bamfile = os.popen("{samtools} view -h {sortedbam}".format(**remove_chroms))
	tempsam = open("{tempsam}".format(**remove_chroms),'w')
	for line in bamfile:
		if(line[0] == "@"):
			tempsam.write(line)
			continue
		elements = line.split('\t')
		if(len(elements[2]) <=5 and elements[2] != "chrM" and elements[2] != "*"):
			tempsam.write(line)
	tempsam.close()
	os.system("{samtools} view -bS {tempsam} > {outputbam}".format(**remove_chroms))
	os.remove("{tempsam}".format(**remove_chroms))
def quality_controls():
	#filter unwanted reads
	remove_junk_chromosome()

	#remove duplicate reads
	
	duplicates = {
	"markduplicate" : param["programs"]["markduplicates"],
	"inputbam" : os.path.join(bamsDir, INPUT + ".sorted.chr.bam"),
	"metricsfile" : os.path.join(bamsDir, INPUT + ".markdup.metrics "),
	"outputbam" : os.path.join(bamsDir, INPUT + ".sorted.chr.nodup.bam"),
	}
	remove_duplicate_CMD = "java -Xmx2g -jar {markduplicate} INPUT={inputbam} \
	METRICS_FILE={metricsfile} OUTPUT={outputbam} REMOVE_DUPLICATES=true \
	ASSUME_SORTED=true".format(**duplicates)
	print(remove_duplicate_CMD)
	os.system(remove_duplicate_CMD)
	
	#mark duplicate reads
	mark_duplicate_CMD = "java -Xmx2g -jar {markduplicate} INPUT={inputbam} \
	METRICS_FILE={metricsfile} OUTPUT={outputbam} REMOVE_DUPLICATES=false \
	ASSUME_SORTED=true".format(**duplicates)
	print(mark_duplicate_CMD)
	os.system(mark_duplicate_CMD)
	
	#remove low quality reads
	filter_param = {
	"samtools" : param["programs"]["samtools"],
	"inputbam" : os.path.join(bamsDir, INPUT + ".sorted.chr.nodup.bam"),
	"outputbam" : os.path.join(bamsDir, INPUT +".sorted.chr.nodup.filt.bam")
	}
	remove_lowqual_CMD ="{samtools}  view -b -h -q 30 {inputbam} > {outputbam}".format(**filter_param)
	print(remove_lowqual_CMD)
	os.system(remove_lowqual_CMD)

	#index bam file
	index_bam_CMD = "{samtools} index {outputbam}".format(**filter_param)
	print(index_bam_CMD)
	os.system(index_bam_CMD)
	
	#qualimap stats
	qualimap_param = {
	"qualimap" : param["programs"]["qualimap"],
	"inputbam" : os.path.join(bamsDir, INPUT +".sorted.chr.nodup.filt.bam"),
	"pdfoutput" : INPUT  + ".filt.sorted.chr.nodup.qualimap.pdf",
	"outputdir" : qualityDir
	}
	qualimap = "{qualimap} bamqc -bam {inputbam} --outdir {outputdir} --outfile {pdfoutput}".format(**qualimap_param)
	
	print(qualimap)
	os.system(qualimap)

def create_insertfile(args):
	bamfile = os.popen("{samtools} view  {inputbam}".format(**args))
	outfile = open("{insert_temp}".format(**args),'w')
	for line in bamfile:
		elements = line.split('\t')
		outfile.write(elements[0]+"\t"+elements[8]+"\n")
	outfile.close()
def shift_cordinates(args):
	infile =  os.popen("paste {tempbed} {insert_temp}".format(**args))
	outfile = open("{shiftedbed}".format(**args),'w')
	for line in infile:
		field = line.strip().split('\t')
		if(field[3].split('/')[0] == field[6] and field[5] == "+"):
			outfile.write(field[0]+"\t"+str(int(field[1].strip())+4)+"\t"+field[2]+\
			"\t"+field[3]+"\t"+field[7]+"\t"+field[5]+"\n")
		elif(field[3].split('/')[0] == field[6] and field[5] == "-"):
			outfile.write(field[0]+"\t"+field[1]+"\t"+str(int(field[2].strip())-5)+\
			"\t"+ field[3]+"\t"+ field[7]+"\t"+field[5]+"\n")
		else:
			print("ERROR "+field[3].split('/')[0],field[6])
	outfile.close()
	os.remove("{tempbed}".format(**args))
	os.remove("{insert_temp}".format(**args))

def peak_calling():
	#shift cordinates
	shift_cordinate_param = {
	"bamtobed" : param["programs"]["bamtobed"],
	"samtools" : param["programs"]["samtools"],
	"inputbam" : os.path.join(bamsDir, INPUT +".sorted.chr.nodup.filt.bam"),
	"tempbed" : INPUT + ".temp.bed",
	"insert_temp" : INPUT +".insert.temp",
	"shiftedbed" : os.path.join(bedsDir, INPUT + ".sorted.chr.nodup.filt.shift.bed")
	}
	temp_bed_CMD = "{bamtobed} -i {inputbam} > {tempbed}".format(**shift_cordinate_param)
	print(temp_bed_CMD)
	os.system(temp_bed_CMD)
	create_insertfile(shift_cordinate_param)
	shift_cordinates(shift_cordinate_param)
	
	#macs peak calling
	macs = {
	"macs_environ" : param["env"]["macs_python_env"],
	"macs" : param["programs"]["macs2"],
	"inputbed" :  os.path.join(bedsDir, INPUT + ".sorted.chr.nodup.filt.shift.bed "),
	"inputname" : INPUT,
	"outputdir" : peaksDir,
	}
	macs_CMD = "export {macs_environ} ; {macs} callpeak --gsize hs -f BED --keep-dup all  \
	--call-summits  --shift -100 --extsize 200 --nomodel --nolambda --verbose 3 \
	--treatment {inputbed} --outdir {outputdir} --name {inputname}".format(**macs)
	print(macs_CMD)
	os.system(macs_CMD)

	#remove blacklisted regions from peaks
	remove_blacklist = {
	"intersectBed" : param["programs"]["intersectBed"],
	"input" : os.path.join(peaksDir, INPUT + "_peaks.narrowPeak"),
	"blacklistfile" : param["files"]["blacklist_regions"],
	"output" : os.path.join(peaksDir, INPUT +"_peaks.narrowPeak.blacklistcleared")
	}
	remove_blacklist_CMD = "{intersectBed} -v -a {input} -b {blacklistfile} > {output}".format(**remove_blacklist)
	print(remove_blacklist_CMD)
	os.system(remove_blacklist_CMD)

def create_bigbed(args):
	infile = open("{inputfile}".format(**args))
	outfile = open("{tempfile}".format(**args),'w')
	for line in infile:
		field = line.strip().split('\t')
		outfile.write(field[0]+"\t"+field[1]+"\t"+field[2]+"\n")
	outfile.close()
	sort_CMD = "sort -k 1,1 -k 2,2n {tempfile} > {sortedpeak} ".format(**args)
	print(sort_CMD)
	os.system(sort_CMD)
	create_bigbed_CMD = "{bedToBigBed} {sortedpeak}  {hg19sizes} {bigbedoutput}".format(**args)
	print(create_bigbed_CMD)
	os.system(create_bigbed_CMD)
	os.remove("{tempfile}".format(**args))
	os.remove("{sortedpeak}".format(**args))

def create_bigwig(args):
	infile = open("{inputbed}".format(**args))
	outfile = open("{tempbed}".format(**args),'w')
	for line in infile:
		fields = line.strip().split("\t")
		if(fields[5] == '-'):
			fields[1] = str(int(fields[2].strip()) - 1)
		end = str(int(fields[1]) + 1)
		outfile.write(fields[0]+"\t"+fields[1]+"\t"+end+"\t"+fields[3]\
		+"\t"+ fields[4]+"\t"+fields[5]+"\n")
	outfile.close()
	sort_CMD = "sort -k 1,1 -k2,2n {tempbed} > {centerbed}".format(**args)
	print(sort_CMD)
	os.system(sort_CMD)
	#extend center bed
	extend_CMD = "{bedtools} slop -i {centerbed} -g {hg19sizes} -b 15 > {15bps_extendedbed}".format(**args)
	print(extend_CMD)
	os.system(extend_CMD)
	os.remove("{tempbed}".format(**args))
	libsize = subprocess.getoutput("{samtools} view -c -F 4 {inputbam}".format(**args))
	scaling_factor = str(1000000 / int(libsize.strip()))
	print(scaling_factor)
	bedgraph_CMD = "{bedtools}  genomecov -i {15bps_extendedbed} -g {hg19sizes} ".format(**args) + \
	"-bg -scale "+ scaling_factor + " > {bedgraphfile}".format(**args)
	print(bedgraph_CMD)
	os.system(bedgraph_CMD)
	bigwig_CMD = "{bedGraphToBigWig} {bedgraphfile} {hg19sizes}  {bigwigoutput}".format(**args)
	print(bigwig_CMD)
	os.system(bigwig_CMD)
	
def visualization():
	create_bigbed_parameters = {
	"inputfile" : os.path.join(peaksDir, INPUT +"_peaks.narrowPeak.blacklistcleared"),
	"tempfile" : os.path.join(tracksDir,  INPUT +".temp"),
	"sortedpeak" : os.path.join(tracksDir,  INPUT +".sortedpeak"),
	"bedToBigBed" : param["programs"]["bedToBigBed"],
	"hg19sizes" : param["files"]["hg19.sizes.txt"],
	"bigbedoutput" : os.path.join(tracksDir,INPUT + ".bb")
	}
	#create bigbed
	create_bigbed(create_bigbed_parameters)
	
	create_bigwig_parameters = {
	"inputbed" : os.path.join(bedsDir, INPUT + ".sorted.chr.nodup.filt.shift.bed"),
	"tempbed" : os.path.join(bedsDir, INPUT + ".temp.center.bed"),
	"centerbed" : os.path.join(bedsDir, INPUT + ".center.bed"),
	"bedtools" : param["programs"]["bedtools"],
	"hg19sizes" : param["files"]["hg19.sizes.txt"],
	"15bps_extendedbed" : os.path.join(bedsDir, INPUT + ".15bpextended.bed") ,
	"samtools" : param["programs"]["samtools"],
	"bedtools" : param["programs"]["bedtools"],
	"inputbam" : os.path.join(bamsDir, INPUT +".sorted.chr.nodup.filt.bam"),
	"bedgraphfile" :  os.path.join(bedsDir, INPUT + ".bedgraph"),
	"bedGraphToBigWig" : param["programs"]["bedGraphToBigWig"],
	"bigwigoutput" : os.path.join(tracksDir,  INPUT +".bw")
	}
	#create bigwig
	create_bigwig(create_bigwig_parameters)

def main():
	check_directory()
	print("\n",time.strftime("%d-%m-%Y %H:%M:%S", time.localtime()),"\n","### START ###",sep='')
	quality_checks()
	mapping()
	quality_controls()
	peak_calling()
	visualization()
	print("\n### DONE ###","\n",time.strftime("%d-%m-%Y %H:%M:%S", time.localtime()),sep='')
if __name__=="__main__":
    main()
