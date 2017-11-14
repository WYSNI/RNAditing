#!/bin/python
#coding:utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import os,subprocess,sys,shutil

#####################################################################################
# Verification of the given percentage

def check_percentage(percentage):
	if percentage >1:
		print("\n ###### Percentage needs to be between 0 and 1 ##########\n")
		sys.exit(2)

#########################################################################################
# Files concatenation

def files_concat(file_out,List):
	concat_file = open(file_out, "w")
	for i in List:
		shutil.copyfileobj(open(i, 'r'), concat_file)
	concat_file.close()

#########################################################################################
# Enumeration of all combination in a List.

def combinations(List_file,n):
    for i in range(len(List_file)):
        v = List_file[i:i+1]
        if n == 1:
            yield v
        else:
            rest = List_file[i+1:]
            for c in combinations(rest, n-1):
                yield v + c

#########################################################################################
# Check path directory

def Exist_dir(path):
	if not os.path.isdir(path):
		print("Error: %s doesn't exist"%path)
		sys.exit(2)

#########################################################################################
# Check path file

def Exist_file(path):
	if not os.path.isfile(path):
		print("Error: %s doesn't exist"%path)
		sys.exit(2)

#####################################################################################
# Common sites beetween two files

def intersection_x2(files,names,path_int,filetype):
	if not os.path.isdir(path_int):
			os.mkdir(path_int)

	if filetype == "Sites":
		
		f1=files[0]
		f2=files[1]
		with open(f1, 'r') as file1:
			with open(f2, 'r') as file2:
				#Common Editing sites beetween file1 and file2 recovery
				same = set(file1).intersection(file2)
	#file name
	name=names[0]+"_"+names[1]
	o_file=name+"_old.txt"
	
	with open(path_int+o_file, 'w') as file_out:
		
		if filetype == "Sites":
			for line in same:
				file_out.write(line)

	#Output file sort
	if filetype == "Sites":
		cmd="sort -k 1,1 -k 2,2n "+path_int+o_file+">"+path_int+name+".bed"

	cmd1="rm "+path_int+o_file
	subprocess.call(cmd,shell =True)
	subprocess.call(cmd1,shell =True)

#####################################################################################
# Common sites beetween three files

def intersection_x3(files,names,path_int,filetype):
	if not os.path.isdir(path_int):
		os.mkdir(path_int)

	if filetype == "Sites":
		
		f1=open(files[0], 'r')
		f2=open(files[1], 'r')
		f3=open(files[2], 'r')
		
		#Common Editing sites beetween file1 and file2 and file3 recovery
		same= set(f1).intersection(f2).intersection(f3)

		f1.close()
		f2.close()
		f3.close()

	# file name
	name=names[0]+"_"+names[1]+"_"+names[2]
	o_file=name+"_old.txt"
	f_o_file=open(path_int+o_file, 'w')

	if filetype == "Sites":
		for line in same:
			f_o_file.write(line)

	f_o_file.close()
	

	#Output file sort
	if filetype == "Sites":
		cmd="sort -k 1,1 -k 2,2n "+path_int+o_file+">"+path_int+name+".bed"
	
	cmd1="rm "+path_int+o_file
	subprocess.call(cmd,shell =True)
	subprocess.call(cmd1,shell =True)

#####################################################################################
# Editing Sites Number

def line_number(file):	
	with open(file, 'r') as file1:
		cpt=0
		for line in file1:
			if line !=" ":
				cpt+=1
	return cpt

#####################################################################################
# Sample or population name recovery

def name_recov (n):
	name=""
	for li in n:
		if name =="":
			name=li
		else:
			name+="_"+li	
	return name

#####################################################################################
# Percentage to integer

def percent(number,total):

	res=float(number)/float(total)*100
	result=round(res, 2)
	return result
