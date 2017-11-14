#!bin/python
# coding=utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import glob,os,subprocess

import functions
import Annotation
import Venn
import Graph
import Clusters

##############################################################################################
# De novo Editing Sites Recovery

def difference(files,path_out):
	if not os.path.isdir(path_out):
		os.mkdir(path_out)
	
	f1=files[0]
	f2=files[1]
	with open(f1, 'r') as file1:
		N=f1[:-4].split("_")[-1]
		with open(f2, 'r') as file2:
			same = set(file1) - set(file2)
			
	same.discard('\n')

	name="Editing_Sites_De_novo_"+N
	o_file=name+"_old.bed"
	
	with open(path_out+o_file, 'w') as file_out:
		cpt=0
		for line in same:
			file_out.write(line)
			cpt+=1

	#Output file sort
	cmd="sort -k 1,1 -k 2,2n "+path_out+o_file+">"+path_out+name+".bed"
	cmd1="rm "+path_out+o_file
	subprocess.call(cmd,shell =True)
	subprocess.call(cmd1,shell =True)

##############################################################################################
# Creating files pairs to compare

def list_difference(list_RNAEditor,list_Res_comp,path_out):
	files=[]
	i=0
	for f in list_RNAEditor:

		files.append(list_RNAEditor[i])
		files.append(list_Res_comp[i])
		difference(files,path_out)
		i+=1
		del files[0:2]

##############################################################################################
# De novo editing sites recovery
# Barplot creation 
# Editing sites Annotation
# Venn diagram creation

def start_analysis(path_RNAEditor, path_Res_comp):
	list_dic=[]
	list_RNAEditor=[]
	list_Res_comp=[]
	print """
	######################################
	### De novo Editing Sites Recovery ###
	######################################
	"""

	l_RNAEditor = glob.glob(path_RNAEditor+'*')
	l_Res_comp = glob.glob(path_Res_comp+'Known/Editing_sites/*')
	
	
	path=path_Res_comp+"De_novo/"
	if not os.path.isdir(path):
		os.mkdir(path)
	
	path_de_novo=path+"Editing_sites/"
	if not os.path.isdir(path_de_novo):
		os.mkdir(path_de_novo)
	
	#Files list creation
	for file1 in l_RNAEditor:
	
		# RNA Editor results files list creation
		if ".bed" in file1 and "Editing_Sites" in file1:

			n=file1[:-4].split("/")[-1].split("_")[-1:]
			name=functions.name_recov(n)
	
			list_RNAEditor.append(file1)

			for file2 in l_Res_comp:

				n=file2[:-4].split("/")[-1].split("_")[-1:]
				NAME=functions.name_recov(n)

				# Comparison results files list creation
				if ".bed" in file2 and NAME == name:
					list_Res_comp.append(file2)
					break

	
	software="Soft_Comp"
	db_name="RNAEditor_REDItools"
	list_difference(list_RNAEditor,list_Res_comp,path_de_novo)

	#Barplot creation
	Graph.barplot(software,path,db_name,list_dic)

	#Editing sites annotation
	Annotation.startAnnotation(path,db_name)

	# Venn Diagram creation
	Venn.start_analysis(path,software)
	
	#Clusters 
	Clusters.Cluster_recovery(path)
	
