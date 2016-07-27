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
#  
#!/usr/bin/python3

import os, sys, json, time

if len(sys.argv)!=3:
	print("Usage: python3 preprocessing.py <jsonfile> <inputfilename>\nExiting!!!!\n")
	exit(0)
if not os.path.exists(sys.argv[1]):
	print("Error: Input jsonfile not found \nExiting!!\n")
	exit(0)

INPUT = sys.argv[2]
def check_directory(param):
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
	### trim adaptors  ###
	trim_adaptors = "trim_galore --paired --nextera "+ \
	param["folders"]["rawreads"]+"/"+ INPUT +"_1.fastq "+ \
	param["folders"]["rawreads"]+"/"+ INPUT +"_2.fastq "+ \
	"-- output_dir "+ param["folders"]["reads"]
	
	### fastqc quality analyzis ###
	fastqc = "fastqc " +\
	param["folders"]["reads"]+ "/" + INPUT + "_1_val_1.fq " +\
	param["folders"]["reads"]+ "/" + INPUT + "_2_val_2.fq -o "+ \
	param["folders"]["quality"]

	print(trim_adaptors)
	print(fastqc)
	#os.system(trim_adaptors)
	#os.system(fastqc)

def main():
	with open("preprocessing_setup.json") as f:
		param = json.load(f)
	check_directory(param)
	print("\n",time.strftime("%d-%m-%Y %H:%M:%S", time.localtime()),"\n","### START ###",sep='')
	
	quality_checks(param)
	print("\n### DONE ###","\n",time.strftime("%d-%m-%Y %H:%M:%S", time.localtime()),sep='')
if __name__=="__main__":
    main()
