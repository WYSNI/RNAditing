#!/bin/python
#coding:utf_8

# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import os,subprocess,re,getopt,sys,glob

import Comparison
import Annotation
import functions
import usage
import Softwares_Results_Comp as Soft

######################################################################################################################################
# Check directory

def commande_check(path):
	print("Make sure Analysis are done before use this program and you have use the merge option for all softwares (REDItools and RNAEditor)\n")
	print("If it's ok, the commande line is :\n  RNAditing -E -m -d /path/to/databaseFile -o /path/to/output_directory/ -p 0.75 -c /path/to/RNAEditor/Results/ ####### \n")
	sys.exit(2)

#############################################################################################################
# Database name recovery

def db_name_recuperation(database):
	db=os.path.basename(database)
	d=re.search('[^.]*',db)
	db_name=d.group(0)
	return(db_name)

#############################################################################################################
# Editing Sites knowns on database recovery

def startREDItools(argv,path_REDItools):

	dbSNP=splice=pipeline=perc=False
	List_return=[]
	nb=0

	for arg in argv:
		if arg =="-d" or arg == "--Database" or arg == "-r" or arg == "--Reference" or arg == "-i" or arg == "--Bam" or arg == "-o" or arg == "--output":
			nb+=1

	# Check for all necessary options
	if nb != 4:

		print("\n*** Error Database, BAM, Results folder and Reference options must be given ***\n")
		usage.REDItools_usage()
		sys.exit(2)


	opts,args= getopt.getopt(argv, "c:p:o:v:s:d:i:r:na",["software_comp","percentage","output","variant","splice","Database","Bam","Reference","analyse","all"])

	# Options check
	for opts,args in opts:
		if opts =='-d' or opts == '--Database':
			functions.Exist_file(args)
			database=args

		elif opts =='-i' or opts == '--Bam':
			functions.Exist_dir(args)
			path_bam=args
		
		elif opts =='-o' or opts == '--output':
			functions.Exist_dir(args)
			path_res=args

		elif opts == '-r' or opts == '--Reference':
			functions.Exist_file(args)
			path_ref=args

		elif opts == "-a" or opts =="--all":
			continue

		elif opts == "-p" or opts =="--percentage":
			percentage=float(args)
			perc=True
			functions.check_percentage(percentage)
			

		elif opts == "-n" or opts =="--analyse":
			continue
		elif opts == '-s' or opts == '--SNP':
			functions.Exist_file(args)
			splice_file=args
			splice=True

		elif opts == '-v' or opts == '--variant':
			functions.Exist_file(args)
			dbSNP_file=args
			dbSNP=True
	
		elif opts == '-c' or opts == '--software_comp':
			path_RNAEditor=args+"Results/Editing_sites/"
			functions.Exist_dir(path_RNAEditor)
			pipeline=True
		else:
			print("%s not found"%opts)
			usage.REDItools_usage()
			sys.exit(2)

	# Database name recovery
	db_name=db_name_recuperation(database)
	List_return.append(db_name)

	path_dir=path_res+"Known_"+db_name+"/"
	if not os.path.isdir(path_dir):
		os.mkdir(path_dir)
	
	List_return.append(path_dir)

	path_dir+="Known/"
	if not os.path.isdir(path_dir):
		os.mkdir(path_dir)

	#if percentage option is used
	if perc:
		List_return.append(percentage)

	#if pipeline option is used
	if pipeline:
		List_return.append(path_RNAEditor)
		List_return.append(path_res)

	# BAM files recovery
	bamfiles=path_bam+"bam_files.txt"
	if not os.path.isfile(bamfiles):
		print("\n*** Error: make sure %s exist ***\n"%bamfiles)
		sys.exit(2)

	with open(bamfiles,"r") as files_bam:
		for file in files_bam:
			
			pop=file.split("/")[-3]
			path_out=path_dir+pop+"/"
			if not os.path.isdir(path_out):
				os.mkdir(path_out)

			print("** Starting analyse for %s **"%file)
			filename=file.split("/")[-1].split(".")[0]
		
			if database[-3:] == ".gz" :
				if dbSNP:
					if splice:
						cmd="python "+ path_REDItools+"reditools/REDItoolKnown.py "+ "-f "+ path_ref+ " -l "+ database+ " -o "+ path_out+ " -i "+ file.replace("\n","")+" -F "+ filename+" -d -p -S -K "+dbSNP_file+" -P "+splice_file
					else:
						cmd="python "+ path_REDItools+"reditools/REDItoolKnown.py "+ "-f "+ path_ref+ " -l "+ database+ " -o "+ path_out+ " -i "+ file.replace("\n","")+" -F "+ filename+" -d -p -S -K "+dbSNP_file
				elif splice:
					cmd="python "+ path_REDItools+"reditools/REDItoolKnown.py "+ "-f "+ path_ref+ " -l "+ database+ " -o "+ path_out+ " -i "+ file.replace("\n","")+" -F "+ filename+" -d -p -S -P "+splice_file
				else:
					cmd="python "+ path_REDItools+"reditools/REDItoolKnown.py "+ "-f "+ path_ref+ " -l "+ database+ " -o "+ path_out+ " -i "+ file.replace("\n","")+" -F "+ filename
			
			
				subprocess.call(cmd,shell=True)
		
			else:
				da='.'.join(database.split('.')[:-1])

				if dbSNP:
					if splice:
						cmd="python "+ path_REDItools+"reditools/REDItoolKnown.py "+ "-f "+ path_ref+ " -l "+ database+ " -o "+ path_out+ " -i "+ file.replace("\n","")+"-X -F "+ filename+" -d -p -S -K "+dbSNP_file+" -P "+splice_file
					else:
						cmd="python "+ path_REDItools+"reditools/REDItoolKnown.py "+ "-f "+ path_ref+ " -l "+ database+ " -o "+ path_out+ " -i "+ file.replace("\n","")+"-X -F "+ filename+" -d -p -S -K "+dbSNP_file
				elif splice:
					cmd="python "+ path_REDItools+"reditools/REDItoolKnown.py "+ "-f "+ path_ref+ " -l "+ database+ " -o "+ path_out+ " -i "+ file.replace("\n","")+"-X -F "+ filename+" -d -p -S -P "+splice_file
				else:		
					cmd="python "+ path_REDItools+"reditools/REDItoolKnown.py "+ "-f "+ path_ref+ " -l "+ database+ " -o "+ path_out+ " -i "+ file.replace("\n","")+"-X -F "+ filename
			
			
				subprocess.call(cmd,shell=True)

				database+=".gz"
				dbSNP_file+=".gz"

		print("** Analyse for %s done **\n"%file )
	return(List_return)

#######################################################################################################################################
# Editing Sites knowns on database recovery
# Editing Sites Annotation

def startAnalysis(argv):
	software="REDItools"
	path_REDItools=os.getenv("RNAditing_PATH_REDItools")

	if len(argv) <1:
		print ("\n*** Options must be given ***\n")
		usage.REDItools_usage()
		sys.exit(2)
	
	# Editing Sites knowns on database recovery
	if sys.argv[2] == "-n" or sys.argv[2] == "--analyse":

		List_return=startREDItools(argv,path_REDItools)
		print("\n*** Results here : %s ***\n"%List_return[1])

	# Editing Sites knowns on database recovery / Editing Sites by population / Annotation
	elif sys.argv[2] == "-a" or sys.argv[2] == "--all":

		n=0
		for arg in argv:
			if arg =="-p" or arg == "--percentage":	
				n+=1

		#Check for all necessary options
		if n != 1:

			print("\n*** Error percentage option must be given ***\n")
			usage.REDItools_usage()
			sys.exit(2)

		List_return=startREDItools(argv,path_REDItools)

		# Editing sites samples comparison by population
		Comparison.startAnalysis(List_return[1],List_return[0],software,List_return[2])

		# Editing sites annotation
		Annotation.startAnnotation(List_return[1]+"Results/",List_return[0])
		
		print("##############################################################################################################################")
		print("\n######## REDItools Results : %s ######## \n"%List_return[1])
		print("##############################################################################################################################")
		
		# If the pipeline is called in full
		if len(List_return) == 5:

			Results_REDItools=List_return[1]+"Results/Editing_sites/"

			# Comparison between RNA Editor and REDItools results
			print(Results_REDItools)
			print(List_return[3])
			print(List_return[4])
			Soft.Start_Comparison(Results_REDItools,List_return[3],List_return[4])
		
				
	# Editing sites by population / Editing sites Annotation
	elif sys.argv[2] == "-m" or sys.argv[2] == "--merge":
		nb=0

		for arg in argv:
			if arg == '-d' or arg =='--Database' or arg == '-o' or arg =='--output' or arg == '-p' or arg =='--percentage':
				nb +=1

		#Check for all necessary options
		if nb != 3:
			print ("\n*** Options database, resuts folder and percentage must be given ***\n")
			usage.REDItools_usage()
			sys.exit(2)

		if sys.argv[3] == '-d' or sys.argv[3]=='--Database' :
			database=sys.argv[4]
			functions.Exist_file(database)
			
			if sys.argv[5] == '-o' or sys.argv[5] =='--output':
				path_dir=sys.argv[6]
				functions.Exist_dir(path_dir)
				
				if sys.argv[7] == '-p' or sys.argv[7]=='--percentage' :
					percentage=float(sys.argv[8])
					functions.check_percentage(percentage)
			
				else:
					print("\n ####### Error please make this command:\n RNAditing -E -m -d /path/to/databaseFile -o /path/to/output_directory/ -p 0.75 ####### \n")
					sys.exit(2)
			
			db_name=db_name_recuperation(database)
			path_out=path_dir+"Known_"+db_name+"/"
			
			# Editing sites samples comparison by population
			Comparison.startAnalysis(path_out,db_name,software,percentage)
		
			# Editing sites annotation
			#Annotation.startAnnotation(path_out+"Results/",db_name)
			if len(argv) >8:
				if sys.argv[9] == '-c' or sys.argv[9]=='--software_comp' :
						if sys.argv[10].split("/")[-2] !="Results":
							commande_check(sys.argv[10])

						Results_REDItools=path_out+"Results/Editing_sites/"
						Final_Results=("/").join(path_out.split("/")[:-2])+"/"
						Soft.Start_Comparison(Results_REDItools,sys.argv[10]+"Editing_sites/",Final_Results)

	else:
		usage.REDItools_usage()
		sys.exit(2)	
	
	





