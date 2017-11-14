#!/bin/python
# coding=utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import os,subprocess,re,glob

import Graph
import Venn
import functions
import Clusters


#####################################################################################
# Percentage to integer conversion

def perc_to_int(percentage,nb_sample):
	
	return int(round(percentage*nb_sample))		

#####################################################################################
# Editing sites position by population

def common_editing_sites_recovery(L_file,path_int,path_dir,percentage):
	pop=L_file[0].split("/")[-4]
	
	#Files concatenation
	file_out=path_int+"concat_file_"+pop+".txt"

	#Editing islands files concatenation of a population
	functions.files_concat(file_out,L_file)

	path_dir+="Results/"
	if not os.path.isdir(path_dir+"Editing_sites/"):
		os.mkdir(path_dir+"Editing_sites/")

	path_ed=path_dir+"Editing_sites/"

	# Minimum percentage of a site for it to be considered representative of the population -> integer
 	nb_sample=perc_to_int(percentage,len(L_file))
	
	#Common Editing Sites Recovery
	cmd =("more "+file_out+"|sort -k1,1 -k2,2n|uniq -c|awk '{if ($1 >= %d) print $2\"\t\"$3\"\t\"$4}'>"%nb_sample)+path_ed+"Editing_Sites_"+pop+".bed" 
	subprocess.call(cmd,shell=True)

###########################################################################################################
# Editing sites positions recovery

def listdirectory_positions(List,path_dir,path_out,software,percentage):
	dic={}
	L=0
	L_name=[];L_file=[];L_islands=[]
	
	dirname=path_out.split("/")[-2]
	
	while L <len(List):
		
		path_dir_out=path_out+List[L]+"/"
		
		if software == "REDItools":
			
			path_AG=path_dir_out+"AG/"
			# File with AG lines recovery
			l = glob.glob(path_AG+'/*')
			
		elif software =="RNAEditor":
			
			# Editing_sites positions recovery
			path_sup=path_dir_out+"SUP/"

			if not os.path.isdir(path_sup):
				os.mkdir(path_sup)

			l = glob.glob(path_dir_out+'/*')
			
			# For file in each sample results directory
			for file in l: 
				end=(".").join(file.split(".")[-2:])

				n=("_").join(file.split("/")[-1].split(".")[:-2]).split("_")[:-1]
				name=functions.name_recov(n)
				
				# RNA Editor editing sites results files recovery
 				if end =="editingSites.gvf":
					
					# Positions recovery
					cmd="grep -v \"#\" "+file+"|awk '{print \"chr\"$4\"\t\"$8\"\t\"$8}' |sort -k 1,1 -k 2,2n |uniq>"+path_sup+name+"_positions.AG.txt"
					subprocess.call(cmd,shell=True)

				if end =="editingIslands.bed":

					cmd="grep -v \"#\" "+file+"|awk '{print \"chr\"$1\"\t\"$2\"\t\"$3}' |sort -k 1,1 -k 2,2n |uniq>"+path_sup+name+"_positions.islands.txt"
					subprocess.call(cmd,shell=True)
					
			l = glob.glob(path_sup+'/*')


		for file in l:

			if "positions.AG.txt" in file:
				# Sample name recovery
				n=(file.split("/")[-1]).split("_")[:-1]
				name=functions.name_recov(n)

				L_file.append(file)
				
				# Dictionary with Editing sites number for each Samples
				dic[name]=functions.line_number(file)

			if "positions.islands.txt" in file:
				L_islands.append(file)
		L+=1

	if not os.path.isdir(path_dir+"Results/"):
		os.mkdir(path_dir+"Results/")
	path_int=path_dir+"Results/"

	if not os.path.isdir(path_int+"Tmp/"):
		os.mkdir(path_int+"Tmp/")
	path_int+="Tmp/"

	# Common Editing Sites Recovery
	Editing_sites_out=common_editing_sites_recovery(L_file,path_int,path_dir,percentage)
	
	os.system("rm -rf "+path_int)

	return dic


###########################################################################################################
# Editing sites recovery	

def listdirectory_AG(List,path_out):
		L=0
		while L <len(List):

			path_dir_out=path_out+List[L]+"/"
			path_AG=path_dir_out+"AG/"
			l = glob.glob(path_AG+'/*')
  
			for file in l: 
				if file.split("_")[-1]=="AG.txt":
					
					n=(file.split("/")[-1]).split("_")[:-1]
					name=functions.name_recov(n)

					#positions recovery
					cmd="cat "+file+"|awk '{print \"chr\"$1\"\t\"$2\"\t\"$2}'|sort -k 1,1 -k 2,2n|uniq >"+path_AG+name+"_positions.AG.txt"
					subprocess.call(cmd,shell=True)
			L+=1


###########################################################################################################		
# AG lines recovery

def AG_recovery(file,name,path_AG):
	
	cmd="cat "+file+"|grep AG |sort -k 1,1 -k 2,2n>"+path_AG+name+"_AG.txt"
	subprocess.call(cmd,shell=True)


###########################################################################################################
# Output table recovery

def listdirectory_outTable(List,path_out):
	L=0
	# Until all samples have been analyzed
	while L <len(List):
		
		path_dir_out=path_out+List[L]+"/"
	
		if not os.path.isdir(path_dir_out+"AG/"):
			os.mkdir(path_dir_out+"AG/")

		path_AG=path_dir_out+"AG/"
		l = glob.glob(path_dir_out+'/*')  
		
		for file in l: 
			if (file.split("/")[-1]).split("_")[0]=="outTable":
			# Sample name recovery
				n=("_").join(file.split("/")[-2].split("_")[:-2]).split("_")[1:]
				name=functions.name_recov(n)

			# Editing sites recovery
				AG_recovery(file,name,path_AG)
		
		L+=1
	
############################################################################################################
# Process list for each software
	
def process_list(List,path_dir,path_out,software,percentage):

		if software == "REDItools":

			listdirectory_outTable(List,path_out)
			listdirectory_AG(List,path_out)
			dic=listdirectory_positions(List,path_dir,path_out,software,percentage)

		elif software =="RNAEditor":

			dic=listdirectory_positions(List,path_dir,path_out,software,percentage)

		return dic

#############################################################################################################
# Per population editing sites recovery
# Barplot creation
# Venn diagram creation

def startAnalysis(path_dir,db_name,software,percentage):
	dic={}
	List_dic=[]
	if software == "REDItools":
		path_out=path_dir+"Known/"

	elif software == "RNAEditor":
		path_out=path_dir+"Results_RNAEditor/"

	if not os.path.isdir(path_out):
		os.mkdir(path_out)

	print("\n*** Files Comparisons ***\n")

	# For each population directory
	for directory in os.listdir(path_out):
		List=[]
		path_d=path_out+directory+"/"

		# for each sample of this population
		for direc in os.listdir(path_d):
			List.append(direc)
		
		# Analysis by population
		print("# %s files Analysis # "%directory)
		dic=process_list(List,path_dir,path_d,software,percentage)
		List_dic.append(dic)

	print("\n*** Barplot and Stat files Creation ***\n")

	path_results=path_dir+"Results/"
	# Barplot Creation
	Graph.barplot(software,path_results,db_name,List_dic)

	# Venn Diagram Creation
	#Venn.start_analysis(path_results,software)	

