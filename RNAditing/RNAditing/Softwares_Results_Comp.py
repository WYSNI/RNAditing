#!/bin/bash
# coding=utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import os, subprocess, re, glob, sys

import functions
import Graph
import Annotation
import Venn
import De_novo
import Clusters


######################################################################################################################################
# Dictionary creation

def Make_dic(path):
	dic={}
	l=glob.glob(path+'/Editing_Sites_*')
	for file in l :
		name=""
		if file [-4:]==".bed":
			
			n =file[:-4].split("/")[-1].split("_")[-1:]
			for li in n:

				if name =="":
					name=li
				else:
					name+="_"+li	
			dic[name]=file
	return dic
######################################################################################################################################
# Comparison between RNAEditor and REDItools Results and De novo Editing sites Recovery

def Start_Comparison(path_REDItools, path_RNAEditor,path_out):
	list_dic=[]

	print """
	##########################################################
	### Comparison between RNAEditor and REDItools Results ###
	##########################################################
	"""

	
	path_results=path_out+"Results/"

	if not os.path.isdir(path_results):
		os.mkdir(path_results)

	path_results+="Known/"
	if not os.path.isdir(path_results):
		os.mkdir(path_results)

	path_results+="Editing_sites/"
	if not os.path.isdir(path_results):
		os.mkdir(path_results)

	# Files Recovery
	dic_RED=Make_dic(path_REDItools)
	dic_RNAEd=Make_dic(path_RNAEditor)
	
	# For each_population
	for key in dic_RED.keys():
		List_files=[]
		List_names=[]
		
		# Same files recovery (/path/to/REDitools/Results/Editing_Sites/Editing_Sites_pop1.bed /path/to/RNAEditor/Results/Editing_Sites/Editing_Sites_pop1.bed)
		List_files.append(dic_RED[key])
		List_files.append(dic_RNAEd[key])
	
		List_names.append("Editing_Sites_REDItools_"+key)
		List_names.append("RNAEditor_"+key)
		
		
		# Comparison between 2 files
		functions.intersection_x2(List_files,List_names,path_results,"Sites")
	
	db_name="RNAEditor_REDItools"
	software="Soft_Comp"

	Known_path=("/").join(path_results.split("/")[:-2])+"/"
	print("*** Graph Creation ***\n")	
	# Barplot creation

	Graph.barplot(software,Known_path,db_name,list_dic)

	# Venn diagram creation
	#Venn.start_analysis(Known_path,software)
	
	# Editing sites annotation
	Annotation.startAnnotation(Known_path,db_name)
	
	Results_path=("/").join(path_results.split("/")[:-3])+"/"

	
	# Clusters recovery
	Clusters.Cluster_recovery(Known_path)

	# De novo Editing sites recovery
	#De_novo.start_analysis(path_RNAEditor,Results_path)

	print("##############################################################################################################################")
	print("\n######## Final Results : %s ######## \n"%Results_path)
	print("##############################################################################################################################")

