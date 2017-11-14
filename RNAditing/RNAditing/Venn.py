#!/bin/python
# coding:utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import glob,os,difflib

import Comparison
import functions
import Graph

##########################################################################################################
# Population specific editing sites recovery 

def diff(file1,file2,name,path):
	f1=open(file1,"r")
	
	f2=open(file2,"r")
	
	#Only different lines beetween file1 and file2 are kept
	diff = difflib.ndiff(f1.readlines(), f2.readlines())
	line = ''.join(x[2:] for x in diff if x.startswith('- '))

	o_file=name+"_only.bed"
	with open(path+o_file, 'w') as file_out:
		file_out.write(line)

	f1.close()
	f2.close()

	# A file name_only.bed is created
	return o_file

##########################################################################################################
# Intersection and Population editing sites recovery for a Venn Diagram creation

def start_analysis(path_dir,software):
	List_files=[]
	List_name=[]

	print("\n# Venn Diagram Creation #\n")
	
	path_ed=path_dir+"Editing_sites/"
	
	l = glob.glob(path_ed+'/*')  
	for file in l: 
		if file[-4:]==".bed":
			List_files.append(file)
			
			# Population name recovery
			if software != "Soft_Comp":
				n=file.split(".")[-2].split("/")[-1].split("_")[2:]
			else:
				n=file[:-4].split("/")[-1].split("_")[-1:]

			name=functions.name_recov(n)
			List_name.append(name)
			
	if not os.path.isdir(path_ed+"Venn/"):
		os.mkdir(path_ed+"Venn/")
	path_int=path_ed+"Venn/"


	#Combination calculation for 3/5 files
	
	enum_files_x3=functions.combinations(List_files,3)
	enum_names_x3=functions.combinations(List_name,3)
	

	for files,names in zip(enum_files_x3,enum_names_x3):

		dic={}
		
		#Nombre de ligne total des 3 fichiers MN GC Na
		dic["name1"]=names[0]
		dic["name2"]=names[1]
		dic["name3"]=names[2]
		dic["name1_number"]=functions.line_number(files[0])
		dic["name2_number"]=functions.line_number(files[1])
		dic["name3_number"]=functions.line_number(files[2])
		
		n1=names[0]+"_"+names[1]+"_"+names[2]
		dic["name3_number"]=functions.line_number(files[2])
		folder=path_int+n1+"/"
	
		
		dic["filename"]=(folder+"Venn_"+n1+".png")

		if not os.path.isdir(folder):
			os.mkdir(folder)
		
		#Intersection beetween 3 files
		functions.intersection_x3(files,names,folder,"Sites")
		
		dic["1"]=functions.line_number(folder+n1+".bed")
	
		#Combination calculation for 2/3 files
		enum_files_x2=functions.combinations(files,2)
		enum_names_x2=functions.combinations(names,2)
		
		i=2
		for file,name in zip(enum_files_x2,enum_names_x2):
			
			#Intersection beetween 2 files
			functions.intersection_x2(file,name,folder,"Sites")

			n2=name[0]+"_"+name[1]

			# Differents sites beetween 2 files
			n2=diff(folder+n2+".bed",folder+n1+".bed",n2,folder)

			dic[str(i)]=functions.line_number(folder+n2)
			i+=1

		#Venn diagram Creation
		Graph.Venn_Creation(dic)




