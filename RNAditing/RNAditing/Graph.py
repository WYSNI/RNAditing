#!/bin/python
# coding:utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import glob,os,re,subprocess
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt

import functions

#############################################################################################################################################
# Barplot and information files creation 

def barplot(software,path_dir,db_name,List_dic):

	List_nbr=List_name=""
	i=1

	print("# Barplot Creation #\n")

	path_ed=path_dir+"Editing_sites/"	

	path_stat=path_ed+"Stats_"+db_name+".csv"
	f_stat=open(path_stat,"w")
	
	# Header of Stats files writing
	f_stat.write("Name"+"\t"+"Start Editing Sites Number"+"\t"+"End Editing Sites Number"+"\t"+"Editing Sites losses"+"\t"+"Loss ratio"+"\n")

	# Editing Sites files recovery
	l = glob.glob(path_ed+'/*.bed')  

	# For each Editing sites files
	for file in l:
		
		# Population name recovery
		if software == "Soft_Comp":
			n=file[:-4].split("/")[-1].split("_")[-1:]
			name=functions.name_recov(n)
		else:
			n=(file.split(".")[-2].split("/")[-1].split("_")[2:])
			name=functions.name_recov(n)

		# Number of editing sites
		Nbr=functions.line_number(file)
		
		# For each dictionary
		for dic in List_dic:

			#For each dictionary key
			for key,value in zip(dic.keys(),dic.values()):
				# Information calculation

				if key.split("_")[0]==name or key.split("_")[1]==name:
					loss=value-Nbr
					rat=float(loss)/float(value)
					ratio=format(rat, '.2f')
					f_stat.write(key+"\t"+str(value)+"\t"+str(Nbr)+"\t"+str(loss)+"\t"+str(ratio)+"\n")
				
		
		# List of number and names in string format for the R Script
		if i==1:
			List_nbr=str(Nbr)
			List_name=str(name)
		else:
			List_nbr=List_nbr+","+str(Nbr)
			List_name=List_name+","+str(name)
		
		i+=1

	f_stat.close()

	filename=path_ed+"Editing_Sites_"+db_name+".png"
	
	#Using a R Script for the barplot creation 
	path_RNAditing=os.getenv("RNAditing_PATH_RNAditing")
	cmd="Rscript "+path_RNAditing+"barplot.R "+List_nbr+" "+List_name+" "+filename+" "+db_name+"2>"+"trash.txt"
	subprocess.call(cmd,shell=True)
	
	cmd="rm trash.txt"
	subprocess.call(cmd,shell=True)
	
	print("2 files created")
	print("- %s"%path_stat)
	print("- %sEditing_Sites_%s.png\n"%(path_ed,db_name))

#############################################################################################################################################
# Venn diagram creation

def Venn_Creation(dic):
	

	# Rest calculation. (For Venn left, right and bottom parts respectively)
	var5=dic.get("name1_number")-dic.get("2")-dic.get("3")-dic.get("1")
	var6=dic.get("name2_number")-dic.get("2")-dic.get("4")-dic.get("1")
	var7=dic.get("name3_number")-dic.get("3")-dic.get("4")-dic.get("1")

	# Total of sites
	total=dic.get("1")+dic.get("2")+dic.get("3")+dic.get("4")+var5+var6+var7

	# Percentages for each part of the Venn Diagram
	p1=functions.percent(dic.get("1"),total)
	p2=functions.percent(dic.get("2"),total)
	p3=functions.percent(dic.get("3"),total)
	p4=functions.percent(dic.get("4"),total)
	p5=functions.percent(var5,total)
	p6=functions.percent(var6,total)
	p7=functions.percent(var7,total)
	
	# Subset sizes
	s = (
	    5,    # Abc 100 left
	    5,    # aBc	010 right
	    5,    # ABc 110 center top
	    5,    # abC 001 bottom
	    5,    # AbC 101 center left
	    5,    # aBC 011 center right
	    1,    # ABC 111 center
	)

	# Venn diagram creation
	fig=plt.figure()
	c = venn3_circles(subsets=s, linestyle='solid')
	v = venn3(subsets=s,set_labels = (dic.get("name1")+"\n("+str(dic.get("name1_number"))+")", dic.get("name2")+"\n("+str(dic.get("name2_number"))+")", dic.get("name3")+"\n("+str(dic.get("name3_number"))+")"))
	
	v.get_label_by_id('111').set_text(str(dic.get("1"))+"\n("+str(p1)+"%)")
	v.get_label_by_id('110').set_text(str(dic.get("2"))+"\n("+str(p2)+"%)")
	v.get_label_by_id('101').set_text(str(dic.get("3"))+"\n("+str(p3)+"%)")
	v.get_label_by_id('011').set_text(str(dic.get("4"))+"\n("+str(p4)+"%)")
	v.get_label_by_id('100').set_text(str(var5)+"\n("+str(p5)+"%)")
	v.get_label_by_id('010').set_text(str(var6)+"\n("+str(p6)+"%)")
	v.get_label_by_id('001').set_text(str(var7)+"\n("+str(p7)+"%)")
	
	
	
	# Venn diagram writing
	
	plt.title("Editing Sites Numbers")
	plt.draw()
	plt.savefig(dic.get("filename"))
	plt.close(fig)

	print("- %s created\n"%dic.get("filename"))


