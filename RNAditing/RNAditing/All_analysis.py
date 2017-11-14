#!/bin/python
#coding:utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import subprocess,os,getopt,sys,glob,re

import RNAEditor
import functions
import usage



def startAnalysis(argv):
	nb=0

	#Analysis option must be the first option called
	for arg in argv:
		if not sys.argv[2] == "-a" or sys.argv[2] == "--analysis":
			print("\n################## Error: analysis option must be given in first ##################\n")
			usage.All_Analysis_usage()
			sys.exit(2)
		
		#Check for all necessary options
		if arg =="-c" or arg == "--configure" or arg == "-f" or arg == "--fastq" or arg == "-s" or arg == "--sample" or arg == "-o" or arg == "--output" or arg == "-a" or arg == "--analysis" or arg =="-p" or arg == "--percentages" or arg =="-e" or arg == "--end" :
			nb+=1

	if nb != 7:

		print("\n################## Error: configuration file, fastq folder, output folder, sample file, percentage, sequencing and analysis options must be given ##################\n")
		usage.All_Analysis_usage()
		sys.exit(2)

	
	
	#Editing sites Prediction
	L_out=RNAEditor.startAnalysis(sys.argv[2:])

	#List of commands to call the rest of the pipeline
	print("\n###### Prediction ok: ######\n")
	print("# Path to RNAEditor results directory: %s #"%L_out[0])
	print("# Path to BAM files: %s #\n"%L_out[2])

	print("##### To continue the analysis please make the following command lines #####\n")
	print("$ source deactivate ")
	print("$ source activate REDItools")
	print("$ RNAditing -E -a -i %s -c %s [options]\n"%(L_out[2],L_out[0]))
	 	





	

