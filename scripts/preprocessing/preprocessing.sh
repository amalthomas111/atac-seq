#PBS -S /bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l  mem=10000MB
#PBS -l vmem=10000MB
#PBS -l pmem=10000MB
#PBS -l walltime=100:00:00
#PBS -q cmb
#!/bin/bash
cd "$PBS_O_WORKDIR"


### folders required#####
#reads
#rawreads
#quality
#bams
#beds
#peaks
#tracks
############################
#### scripts required#######
#getinsertsize.R
#ATAC_BAM_shifter_gappedAlign.pl
#############################

printf "start: \n$(date)\n"
## convert paired end sra to fastq ##
#fastq-dump --split-3 rawreads/$INPUT_FILE".sra" -O rawreads
#printf "fastqdump over \n $(date)\n"tqdump over \n $(date)\n"

## Trim adapters ##
printf "Trim galore \n $(date)\n"
trim_galore --paired --quality 20 --nextera rawreads/$INPUT_FILE"_1.fastq" rawreads/$INPUT_FILE"_2.fastq" --output_dir  reads

##Run fastqc ##
printf "Trim_galore over\n$(date)\nChecking quality"
fastqc reads/$INPUT_FILE"_1_val_1.fq" reads/$INPUT_FILE"_2_val_2.fq" -o quality

printf "\n$(date)\n File renaming\n"
mv reads/$INPUT_FILE"_1_val_1.fq" reads/$INPUT_FILE"_1_trimmed.fastq"
mv reads/$INPUT_FILE"_2_val_2.fq" reads/$INPUT_FILE"_2_trimmed.fastq"

###### MAPPING & FILTERING  START ######

## map reads using bowtie2 ##
printf "\nRenaming done\n$(date)\nBowtie \n"
bowtie2 -p 8 -X 2000 --fr --no-discordant --no-mixed --minins 38 \
--met-file bams/$INPUT_FILE".alignmetrics.txt" \
-x "/home/rcf-40/amalthom/staging_work/2.genome/bowti2_index/hg19/hg19" \
-1 reads/$INPUT_FILE"_1_trimmed.fastq" -2 reads/$INPUT_FILE"_2_trimmed.fastq" \
-S bams/$INPUT_FILE".sam"

## sort sam ##
printf "\nBowtie done\n$(date)\nSorting by samtools\n$(date)\n"
samtools view -bS bams/$INPUT_FILE".sam" |samtools sort - bams/$INPUT_FILE".sorted"

## Remove unwanted chromosome ##
printf "\nSorting done \n$(date)\nRemove unwanted chroms \n$(date)\n"
samtools view -h bams/$INPUT_FILE".sorted.bam"\
|awk 'substr($0,1,1) == "@"|| (length($3)<=5 && $3!="chrM" && $3 != "*"){{print $0}}'\
|samtools view -bS - > bams/$INPUT_FILE".sorted.chr.bam"

## Remove duplicates ##
printf "\nRemove unwanted chroms done\n$(date)\nRemove duplicates\n"
java -Xmx2g -jar /home/rcf-40/amalthom/panases_soft/picard_1.122/MarkDuplicates.jar \
INPUT=bams/$INPUT_FILE".sorted.chr.bam" \
METRICS_FILE=bams/$INPUT_FILE".markdup.metrics" \
OUTPUT=bams/$INPUT_FILE".sorted.chr.nodup.bam" \
REMOVE_DUPLICATES=true \
ASSUME_SORTED=true

## Remove low quality reads ##
printf "\nRemove duplicates done\n$(date)\nFilter_bam using samtools\n"
samtools view -b -h -q 30 bams/$INPUT_FILE".sorted.chr.nodup.bam" \
 > bams/$INPUT_FILE".sorted.chr.nodup.filt.bam"

##### MAPPING & FILTERING END #######

##### QUALITY CHECKS START ######

## Insert size and alignment stats ##
printf "\nInsert size\n"
samtools view -f66 bams/$INPUT_FILE".sorted.chr.nodup.filt.bam" |cut -f 9|sed 's/^-//' > $INPUT_FILE"_insertsize.txt"
Rscript getinsertsize.R $INPUT_FILE"_insertsize.txt"
mv $INPUT_FILE"_insertsize.txt"  "hist_"$INPUT_FILE"_insertsize.txt.png" $INPUT_FILE"_insertsize.txt_density.png"quality/

## mapped stats ##
samtools index bams/$INPUT_FILE".sorted.chr.nodup.filt.bam"  bams/$INPUT_FILE".sorted.chr.nodup.filt.bam.bai"
samtools flagstat  bams/$INPUT_FILE".sorted.chr.nodup.filt.bam"  \
	quality/$INPUT_FILE".sorted.chr.nodup.filt.bam.flagstat.qc"

## qualimap ###
$qualimap bamqc -bam bams/$INPUT_FILE".sorted.chr.nodup.filt.bam"  -outdir quality -outfile \
 $INPUT_FILE".sorted.chr.nodup.filt.qualimap.pdf"
## qualimap ##

###### QUALITY CHECKS END ######

## Shifting coordinates to center of cut site ##
printf "\nCentering cordinates\n"
perl ATAC_BAM_shifter_gappedAlign.pl bams/$INPUT_FILE".sorted.chr.nodup.filt.bam" \
bams/$INPUT_FILE".sorted.chr.nodup.filt.shift.temp"
samtools sort bams/$INPUT_FILE".sorted.chr.nodup.filt.shift.temp.bam" \
 bams/$INPUT_FILE".sorted.chr.nodup.filt.shift"
rm bams/$INPUT_FILE".sorted.chr.nodup.filt.shift.temp.bam"
printf "\nCentering done\n$(date)\n"

## cut site center bed file ##
bedtools bamtobed  -i bams/$INPUT_FILE".sorted.chr.nodup.filt.shift.bam" | \
awk 'BEGIN{OFS="\t"}{if($6 == "-") $2=$3-1; print $1, $2, $2+1, $4, $5, $6}'| sort -k 1,1 -k2,2n >  \
beds/$INPUT_FILE".sorted.chr.nodup.filt.shift.center.bed"
bedToBam -i beds/$INPUT_FILE".sorted.chr.nodup.filt.shift.center.bed"  \
-g /home/rcf-40/amalthom/panases_soft/hg19.sizes | samtools sort -  \
bams/$INPUT_FILE".sorted.chr.nodup.filt.shift.center"
samtools index bams/$INPUT_FILE".sorted.chr.nodup.filt.shift.center.bam"

## homer tag directory ##
makeTagDirectory quality/"homertag_"$INPUT_FILE -genome hg19 -keepAll -checkGC  \
 -fragLength 1 bams/$INPUT_FILE".sorted.chr.nodup.filt.shift.center.bam"

## bed file for accessible regions ##
bedtools slop -i beds/$INPUT_FILE".sorted.chr.nodup.filt.shift.center.bed" \
-g /home/rcf-40/amalthom/panases_soft/hg19.sizes \
-b 15 > beds/$INPUT_FILE".15bpshifted.bed" 

bedToBam -i beds/$INPUT_FILE".15bpshifted.bed"  \
-g /home/rcf-40/amalthom/panases_soft/hg19.sizes | samtools sort -  \
bams/$INPUT_FILE".15bpshifted"
samtools index bams/$INPUT_FILE".15bpshifted.bam"

####### PEAK CALLING START ########
## macs ##
export PATH=/home/rcf-40/amalthom/panases_soft/anaconda3/envs/py2.7/bin/:$PATH
macs2 callpeak --gsize hs \
         --treatment beds/$INPUT_FILE".sorted.chr.nodup.filt.shift.bam"  \
         --outdir peaks/ \
         --name $INPUT_FILE \
         --keep-dup all \
         --call-summits \
         --shift -100 --extsize 200 \
         --nomodel --nolambda \
         --verbose 3

##### PEAK CALLING END #########

##### TRACKS START#############

## make bed graph ##

printf "\nBam_to_sorted_bed done\n$(date)\Make_bedgraph_and_bigwig\n"
 librarySize=$(samtools view -c -F 4 bams/$INPUT_FILE".sorted.chr.nodup.filt.shift.bam")
touch quality/$INPUT_FILE"_lib.size."$librarySize
 expr="1000000 / $librarySize"
 scaling_factor=$(echo $expr | bc -l)
 echo "[librarySize: $librarySize]" 
 echo "[scaling_factor: $scaling_factor]"

bedtools genomecov -i beds/$INPUT_FILE".15bpshifted.bed"  \
-g /home/rcf-40/amalthom/panases_soft/hg19.sizes  \
 -bg -scale $scaling_factor > beds/$INPUT_FILE".bedgraph"

## make bigwig ##
 ~/panases_soft/bedGraphToBigWig beds/$INPUT_FILE".bedgraph" \
/home/rcf-40/amalthom/panases_soft/hg19.sizes tracks/$INPUT_FILE".bw"


printf "\nCall macs peak done\n$(date)\nmake peak_track\n"
awk '{{OFS="\t"; print $1, $2, $3}}' peaks/$INPUT_FILE"_peaks.narrowPeak" \
| sort -k 1,1 -k 2,2n > tracks/$INPUT_FILE".temp"

## make bigbed ##
/home/rcf-40/amalthom/panases_soft/bedToBigBed tracks/$INPUT_FILE".temp" \
 /home/rcf-40/amalthom/panases_soft/hg19.sizes tracks/$INPUT_FILE".bb"
rm tracks/$INPUT_FILE".temp"

printf "\nMake peak track done\n$(date)\n"

###### TRACKS END ########
