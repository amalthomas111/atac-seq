# atac-seq pipeline
This repository contains codes for atac-seq processing.
# Pre-processing part
## Trim Adaptors

Using [Trim Galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

Command:
For paired end reads:

'''sh
trim_galore --paired   <read1.fastq> <read2.fastq> --output_dir  <output directory>
'''

Optional arguments for trim_galore):
To remove low quality reads: --quality <cutoff>
## Quality control

In general, the quality of high throughput sequencing data need to be assessed before using them to perform other analyses.
The software [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is designed to do some quality control checks
on raw sequence data coming from high throughput sequencing experiments, by which you can get a quick impression of whether
your data have any obvious problem. The typical command for employing the tool is 
```sh
fastqc <read1_trimmed.fastq> <read2_trimmed.fastq>
```
## Mapping
Using [bowtie2]
```sh
bowtie2 -p 8 -X 2000 --fr --no-discordant --no-mixed --minins 38  
--met-file <alignmetrics.txt> 
-x <bowtie-index>
-1 <read1_trimmed.fastq> -2 <read2_trimmed.fastq> 
-S <output.sam>
```
