#!/bin/python
#coding=utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import sys,os,subprocess,gzip

import functions

#####################################################################################################################################
# Database file creation

#Analysis option must be the first option called
if len(sys.argv) <3:
	print ("\n*** Error:Input AND output files must be given ***")
	print ("python TAB_modification path/to/input_file path/to/output_file\n")
	sys.exit(2)

# Database file recovery
path_in=sys.argv[1]
functions.Exist_file(path_in)

# File out recovery
file_out=sys.argv[2]

# File out path verification 
path_out=("/").join(file_out.split("/")[:-1])
functions.Exist_dir(path_out)
path_tmp=path_out+"tmp.tab"


f_out=open(path_tmp,'w')
if ".gz" in (path_in):
	f_in=gzip.open(path_in,'r')
else:
	f_in=open(path_in,'r')

# Database name recovery
db_name=path_in.split("/")[-1].split("_")[0].split(".")[0]

# For each line in database line
for line in f_in:
	if "chromosome" not in line:
		
		col1=line.split("\t")[0]

		if len(col1)==4:
		
			if db_name=="REDIportal":
				f_out.write(str(col1[-1])+"\t"+str(line.split("\t")[1])+"\t"+str(line.split("\t")[4])+"\n")
			elif db_name=="Human":
				f_out.write(str(col1[-1])+"\t"+str(line.split("\t")[1])+"\t"+str(line.split("\t")[3])+"\n")
			

		elif len(col1)==5:
			if db_name=="REDIportal":
				f_out.write(str(col1[-2:])+"\t"+str(line.split("\t")[1])+"\t"+str(line.split("\t")[4])+"\n")
			elif db_name=="Human":
				f_out.write(str(col1[-2:])+"\t"+str(line.split("\t")[1])+"\t"+str(line.split("\t")[3])+"\n")

		else:	
			if db_name=="REDIportal":	
				f_out.write(str(line.split("\t")[0])+"\t"+str(line.split("\t")[1])+"\t"+str(line.split("\t")[4])+"\n")

			elif db_name=="Human":
				f_out.write(str(line.split("\t")[0])+"\t"+str(line.split("\t")[1])+"\t"+str(line.split("\t")[3])+"\n")
	
f_out.close()

#Sorting positions
cmd = "cat "+path_tmp+ "|sort -k1,1 -k2,2n>"+file_out
cmd1 = "rm "+path_tmp
subprocess.call(cmd,shell=True)
subprocess.call(cmd1,shell=True)
			
		
