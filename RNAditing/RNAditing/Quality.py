#!/bin/python
#coding:utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import glob,os,subprocess,sys

import functions


########################################################################################################################
# Sample name recovery on Sample file
	
def infos(information,file):
	file_in=False
	with open(information,"r") as file1:
			for line in file1:
				
				if file in line:
					file_in=True
					return(line.split("\t")[2].rstrip('\n'))

	if not file_in:
		print("\n ################## This file %s isn't on the sample file %s. Please check this file ##################\n"%(file,information))
		sys.exit(2)
			

########################################################################################################################
# FASTQC TRIMMING FASTQC

def start_analysis(path_fastq,information,end):

	List_fastq=[]
	List_ext=[]
	
	if os.path.isdir(path_fastq+"Quality/"):
		print("\n ################## Quality Controle already made ##################\n")
		return path_fastq+"Quality/TRIM/"
	
	# Hierarchy creation
	if not os.path.isdir(path_fastq+"Quality/"):
		os.mkdir(path_fastq+"Quality/")
	path_quality=path_fastq+"Quality/"
	if not os.path.isdir(path_quality+"FASTQC_1/"):
		os.mkdir(path_quality+"FASTQC_1/")
	path_fastqc1=path_quality+"FASTQC_1/"

	if not os.path.isdir(path_quality+"TRIM/"):
		os.mkdir(path_quality+"TRIM/")
	path_trim=path_quality+"TRIM/"
	
	if not os.path.isdir(path_quality+"FASTQC_2/"):
		os.mkdir(path_quality+"FASTQC_2/")
	path_fastqc2=path_quality+"FASTQC_2/"


	l=glob.glob(path_fastq+"/*/*")
	pre_fastq=""
	
	print("\n### Quality control ###\n")

	# For each file
	for file in l:
		
		
		# List of FASTQ creation
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

			#FASTQ name recovery
			fq=file.split("_")[:-1]
			fastq=functions.name_recov(fq)
			if pre_fastq != fastq:
				List_fastq.append(fastq)
			pre_fastq=fastq
			

	#For each FASTQ
	for file,ext in zip(List_fastq,List_ext):
		pop = file.split("/")[-2]
		if not os.path.isdir(path_fastqc1+pop+"/"):
			os.mkdir(path_fastqc1+pop+"/")

		path_pop=path_fastqc1+pop+"/"

		name=infos(information,("_").join(file.split("_")[:-1]))
		
		############# FASTQC #############
	
		name1=name+"_1"
		f1=file+"_1"
		f1_ext=f1+ext
		f_name1=f1.split("/")[-1]
		
		print("\n### Fastqc for %s ###\n"%f1_ext)

		cmd="fastqc "+f1_ext+" -o "+path_pop
		subprocess.call(cmd,shell = True)
		
		cmd2= "mv "+path_pop+f_name1+"_fastqc.html "+path_pop+name1+"_fastqc.html"
		subprocess.call(cmd2, shell=True)

		#if paired end
		if end =="paired":
			name2=name+"_2"
			f2=file+"_2"
			f2_ext=f2+ext
			f_name2=f2.split("/")[-1]
			print("\nFastqc for %s\n"%f2_ext)

			cmd1="fastqc "+f2_ext+" -o "+path_pop
			subprocess.call(cmd1,shell = True)

			cmd3= "mv "+path_pop+f_name2+"_fastqc.html "+path_pop+name2+"_fastqc.html"
			subprocess.call(cmd3, shell=True)

		############## TRIMMOMATIC ##############
		if not os.path.isdir(path_trim+pop+"/"):
			os.mkdir(path_trim+pop+"/")

		path_pop_trim=path_trim+pop+"/"
		res= os.popen("which trimmomatic")
		path_trimmomatic=res.readlines()[0]
		path_trimmomatic=("/").join(path_trimmomatic.split("/")[:-2])+"/share/trimmomatic/adapters/"
		
		print("Trimming for %s.fastq"%file)

		if end =="single":
			cmd4= "trimmomatic SE -phred33 "+f1_ext+" "+path_pop_trim+name1+ext+" ILLUMINACLIP:"+path_trimmomatic+"TruSeq2-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

		else:
			cmd4= "trimmomatic PE -phred33 "+f1_ext+" "+f2_ext+" "+path_pop_trim+name1+ext+" "+path_pop_trim+name1+".unpaired.fq "+ path_pop_trim+name2+ext+" "+path_pop_trim+name2+".unpaired.fq ILLUMINACLIP:"+path_trimmomatic+"TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

			cmd5= "rm "+path_pop_trim+name1+".unpaired.fq "
			cmd6= "rm "+path_pop_trim+name2+".unpaired.fq "
			
		subprocess.call(cmd4,shell=True)

		#Erase unpaired files

		if end =="paired":
			subprocess.call(cmd5,shell=True)
			subprocess.call(cmd6,shell=True)

		############# FASTQC #############

		if not os.path.isdir(path_fastqc2+pop+"/"):
			os.mkdir(path_fastqc2+pop+"/")

		path_pop=path_fastqc2+pop+"/"
		f1=path_pop_trim+name1
		f1_ext=f1+ext

		print("\n### Fastqc for %s ###\n"%f1_ext)

		cmd7="fastqc "+f1_ext+" -o "+path_pop
		subprocess.call(cmd7,shell=True)

		# if paired end
		if end =="paired":
			f2=path_pop_trim+name2
			f2_ext=f2+ext
			print("\n### Fastqc for %s ###\n"%f2_ext)
			
			cmd8="fastqc "+f2_ext+" -o "+path_pop
			subprocess.call(cmd8,shell=True)


	return path_trim
		

########################################################################################################################		
				

if __name__ =="__main__":

	path_fastq="/media/bioinf/faef2f7b-527b-4b8a-859e-72ad73bd98bc/TEST/test/"
	information="/media/bioinf/faef2f7b-527b-4b8a-859e-72ad73bd98bc/TEST/fichier_test.txt"
	end="paired"
	start_analysis(path_fastq, information, end)
