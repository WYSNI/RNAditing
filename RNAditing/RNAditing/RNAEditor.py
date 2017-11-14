#!/bin/python
#coding:utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

#Make sure paired end files have the same name nameOfFastq_1.fastq nameOfFastq_2.fastq

import subprocess,os,getopt,sys,glob


import Comparison
try:
	import Quality
	import Annotation
	import functions
	import usage
except:
	print("\n#### Error: Open the good Conda environnement #####")
	print("#### source activate RNAEditor or REDItools #####\n")
	sys.exit(2)

#########################################################################################
# Hierarchy modification of RNA Editor results

def hierarchy_modification(path_fastq):
	path_out=("/").join(path_fastq.split("/")[:-4])+"/RNAEditor/"
	if not os.path.isdir(path_out):
		os.mkdir(path_out)

	path_dir=path_out+"Results_RNAEditor/"
	if not os.path.isdir(path_dir):
		os.mkdir(path_dir)

	l= glob.glob(path_fastq+"*/rnaEditor/*")

	for folder in l :
		pop=folder.split("/")[-3]
		path_pop=path_dir+pop+"/"
		if not os.path.isdir(path_pop):
			os.mkdir(path_pop)
		cmd= "mv "+folder+" "+path_pop
		subprocess.call(cmd,shell=True)
		
	return path_out

#########################################################################################
# Quality controle and Editing sites prediction  


def startRNAEditor(argv,path_RNAEditor,pipeline):
	List_ext=[]
	List_fastq=[]
	L_out=[]

	merge=False
	nb=0

	
	for arg in argv:

		if arg =="-c" or arg == "--configure" or arg == "-f" or arg == "--fastq" or arg == "-s" or arg == "--sample" or arg == "-e" or arg == "--end":
			nb+=1

	#Check if all option needed are present
	if nb != 4:

		print("\n ###### Error configuration file, fastq folder, sequencing and samples options must be given ##########\n")
		usage.RNAEditor_usage()
		sys.exit(2)

	try:
		opts,args=getopt.getopt(argv,"e:p:s:o:c:f:aimn",["end","percentage","sample","output","configure","fastq","all","GUI","merge","analyse"])
	except getopt.GetoptError as err:
		print str(err)
		usage.RNAEditor_usage()
		sys.exit(2)

	# Options check
	for opts,args in opts:
		if opts == "-c" or opts =="--configure":
			functions.Exist_file(args)
			path_config_file=args

		elif opts == "-f" or opts =="--fastq":
			functions.Exist_dir(args)
			path_fastq=args

		elif opts == "-a" or opts =="--all":
			merge=True
			continue
		elif opts == "-n" or opts =="--analyse":
			continue

		elif opts == "-o" or opts =="--output":
			functions.Exist_dir(args)
			bam_path=args
			pipeline=True
	
		elif opts == "-p" or opts =="--percentage":
			percentage=float(args)
			functions.check_percentage(percentage)

		elif opts == "-s" or opts =="--sample":
			functions.Exist_file(args)
			information=args
			
		elif opts == "-e" or opts =="--end":
			if args == "single" or args == "paired":
				sequencing=args
			else:
				print("\n#### Error: -e argument must be paired or single ####\n")
				functions.RNAEditor_usage()
				sys.exit(2)

		else:
			usage.RNAEditor_usage()
			sys.exit(2)

	###### QUALITY CONTROLE #####
	path_fastq=Quality.start_analysis(path_fastq,information,sequencing)
	
	if pipeline:
		if not os.path.isdir(bam_path+"BAM/"):
			os.mkdir(bam_path+"BAM/")
		bam_path+="BAM/"

	l=glob.glob(path_fastq+"/*/*")
	
	pre_fastq=""
	for file in l :
		# FASTQ files recovery
		if ".fastq" in file or ".fq" in file:
			
			#Extension recovery
			ext=file.split("/")[-1].split(".")[1:]
			extension=""
			for li in ext:
				if extension =="":
					extension="."+li
				else:
					extension+="."+li

			List_ext.append(extension)
			#fastq prefix
			fq=file.split("_")[:-1]
			fastq=functions.name_recov(fq)

			#List of prefix
			if pre_fastq != fastq:
				List_fastq.append(fastq)
			pre_fastq=fastq
		

	# For each file in list
	for file,ext in zip(List_fastq,List_ext):
		fastq_1=file+"_1"+ext
		fastq_2=file+"_2"+ext
	
		#if paired end
		if sequencing == "paired":

			cmd="python " +path_RNAEditor+"RNAEditor.py"+ " -c " +path_config_file+ " -i " +fastq_1+" "+fastq_2 
			
		else:	
			
			cmd="python " +path_RNAEditor+"RNAEditor.py"+ " -c " +path_config_file+ " -i " +fastq_1 

		
		subprocess.call(cmd,shell=True)

		path_pop=("/").join(fastq_1.split("/")[:-1])
		pop=path_pop.split("/")[-1]
		sample=fastq_1.split(".")[0].split("/")[-1]
		path_file=path_pop+"/rnaEditor/Res_"+sample+"/"
		"""
		if pipeline:
			
			if not os.path.isdir(bam_path+pop+"/"):
				os.mkdir(bam_path+pop+"/")
			bam_pop=bam_path+pop+"/"
			cmd= "cp "+path_file+sample+".noDup.realigned.recalibrated.bam "+bam_pop
			subprocess.call(cmd, shell=True)
		"""	
	#Hierarchy_modification
	path_out=hierarchy_modification(path_fastq)
	L_out.append(path_out)

	#BAM files on bam_files.txt	
	files=glob.glob(path_out+"Results_RNAEditor/*/*/*.noDup.realigned.recalibrated.bam")
	

	with open (bam_path+"bam_files.txt","w") as bam_files:
		for file in files: 
			bam_files.write(file+"\n")
		

	# If prediction is made in full
	if merge:
		L_out.append(percentage)

	# If the pipeline is totally called
	if pipeline:
		L_out.append(bam_path)
	
	return L_out

	
#########################################################################################
# Quality controle,Editing sites prediction
# Editing sites by population
# Annotation  

def startAnalysis(argv):

	L_out=[]
	pipeline=False
	path_RNAEditor=os.getenv("RNAditing_PATH_RNAEditor")

	software= "RNAEditor"
	db_name= "RNAEditor"

	if len(argv) <1:
		usage.RNAEditor_usage()
		sys.exit(2)

	#Interface for Prediction
	if sys.argv[2] == "-i" or sys.argv[2] == "--GUI":

		os.chdir(path_RNAEditor)
		if not os.path.isfile("configuration.txt"):
			print("Error configuration file doesn't exist. Make sure configuration.txt is in the same directory as the software")
			sys.exit()

		cmd=" python " +path_RNAEditor+"RNAEditor.py"
		subprocess.call(cmd,shell=True)
		sys.exit(1)

	# Prediction
	elif sys.argv[2] == "-n" or sys.argv[2] == "--analyse":

		L_out=startRNAEditor(argv,path_RNAEditor, pipeline)
		print( "Results are here : %s"%L_out[0])
	
	# Prediction and Editing sites by population
	elif sys.argv[2] == "-a" or sys.argv[2] == "--all":
		n=0
		for arg in argv:
			
			if arg == "-p" or arg == "--percentage":
				n+=1
	
		#Check if all option needed are present
		if n != 1:

			print("\n*** Error percentage option must be given ***\n")
			usage.RNAEditor_usage()
			sys.exit(2)

		L_out=startRNAEditor(argv,path_RNAEditor, pipeline)
		print(L_out)

		#Editing sites by population
		Comparison.startAnalysis(L_out[0],db_name,software,L_out[1])

		# Editing sites annotation
		Annotation.startAnnotation(L_out[0],db_name)
		if len(L_out) == 3:
			return L_out
				
	# Editing sites by population
	elif sys.argv[2] == "-m" or sys.argv[2] == "--merge":
		print("\n####### IMPORTANT: Make sure you have previously made python RNAEditor -n or -a before use this command ################\n")
		n =0
		for arg in argv:	
			if arg == "-f" or arg == "--fastq" or arg =="-p" or arg == "--percentage":
				n +=1

		#Check if all option needed are present
		if n != 2:
			print("Error percentage and fastq folder options must be given")
			usage.RNAEditor_usage()
			sys.exit(2)

		if sys.argv[3] == "-f" or sys.argv[3] == "--fastq":
			
			path_fastq=("/").join(sys.argv[4].split("/")[:-2])+"/RNAEditor/"
			functions.Exist_dir(path_fastq)
			
			if sys.argv[5] == "-p" or sys.argv[5] == "--percentage":
				
				percentage=float(sys.argv[6])
				functions.check_percentage(percentage)	
			else:
				print("\n ####### Error please make this command:\n python RNAditing -R -m -f /path/fastq/ -p 0.75 ####### \n")
				sys.exit(2)

		#Editing sites by population
		Comparison.startAnalysis(path_fastq,db_name,software,percentage)

		# Editing sites annotation
		Annotation.startAnnotation(path_fastq,db_name)

	else:
		usage.RNAEditor_usage()
		sys.exit(2)	



	




	
