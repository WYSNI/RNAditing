#!/bin/python
# coding:utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import subprocess,glob,os,re,csv

import functions
import Genes_recovery

###################################################################################################################
# RNA Editing region percentages recovery

def region_recovery(path_annotation,db_name):

	f_out=open(path_annotation+"Editing_sites_region_"+db_name+".summary.csv","w")
	f_out.write("Name\tIntron\tIntergenic\t3\' UTR\t5\' UTR\tNon-coding\tExon\tOther\n")
	
	for file in os.listdir(path_annotation):
		if file.split("_")[-1]=="Annotation.csv":
			nb_intron=nb_UTR3=nb_UTR5=nb_inter=nb_non_coding=nb_other=nb_exon=0
			p_intron=p_UTR3=p_UTR5=p_inter=p_non_coding=p_other=0.0

			# Population name recovery
			n=("_").join(file.split("_")[-2:]).split("_")[:-1]
			name=functions.name_recov(n)
			with open(path_annotation+file,"r") as file_in:
				for line in file_in:
					if line.split("\t")[2]!="Start":
						
						region=line.split("\t")[7]
						if re.match(r"intron.*",region):
							nb_intron+=1
						elif re.match(r"Intergenic.*",region):
							nb_inter+=1
						elif re.match(r"3'\sUTR.*",region):
							nb_UTR3+=1
						elif re.match(r"5'\sUTR.*",region):
							nb_UTR5+=1
						elif re.match(r"non-coding.*",region):
							nb_non_coding+=1
						elif re.match(r"exon.*",region):
							nb_exon+=1
						else :
							nb_other+=1
			Total=nb_intron+nb_inter+nb_UTR3+nb_UTR5+nb_non_coding+nb_exon+nb_other

			# Percentage Calculation
			p_intron=functions.percent(nb_intron,Total)
			p_UTR3=functions.percent(nb_UTR3,Total)
			p_UTR5=functions.percent(nb_UTR5,Total)
			p_inter=functions.percent(nb_inter,Total)
			p_non_coding=functions.percent(nb_non_coding,Total)
			p_exon=functions.percent(nb_exon,Total)
			p_other=functions.percent(nb_other,Total)
			
			#Writing to a file
			f_out.write(name+"\t"+(str(p_intron)).replace(".",",")+"\t"+(str(p_inter)).replace(".",",")+"\t"+(str(p_UTR3)).replace(".",",")+"\t"+(str(p_UTR5)).replace(".",",")+"\t"+(str(p_non_coding)).replace(".",",")+"\t"+(str(p_exon)).replace(".",",")+"\t"+(str(p_other)).replace(".",",")+"\n")

			
	f_out.close()
	print ("\t- %sEditing_sites_region_%s.summary.csv created\n"%(path_annotation,db_name))


########################################################################################################################################
# Editing sites annotation, edited genes and region percentage recovery
			
def startAnnotation(path,db_name):

	print("*** Annotation ***\n")

	path+="Editing_sites/"
	path_annotation=path+"Annotation/"

	if not os.path.isdir(path_annotation):
		os.mkdir(path_annotation)
	
	
	#For each Editing sites file
	for file in os.listdir(path):
		
		if file[-4:]==".bed":

			n=file[:-4].split("_")[-1:]
			name=functions.name_recov(n)
			
			# Annotation

			file_out=path_annotation+"Editing_sites_"+name+"_Annotation.csv"
			cmd ="annotatePeaks.pl "+path+file+ " hg19 -annStats "+path_annotation+"Annotation_"+name+"_Stats.csv >"+file_out+" 2>"+path_annotation+"trash.txt"
			subprocess.call(cmd,shell=True)

			print("\t- Annotation for %s file OK: %s and %sAnnotation_%s_Stats.csv created\n\n"%(name,file_out,path_annotation,name))
		
			# Edited genes recovery
			Genes_recovery.conversion(file_out,name,path_annotation)
	
	
	cmd="rm "+path_annotation+"trash.txt"
	subprocess.call(cmd, shell=True)
	
	print("*** RNA Editing Region percentages recovery ***\n")
	region_recovery(path_annotation,db_name)

