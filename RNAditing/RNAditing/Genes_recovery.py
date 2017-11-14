#!/bin/python
#coding:utf_8

# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

import subprocess,mygene,os


def ncbi_to_refgene(line):
	refgene_ok=False
	list_info=[]
	mg = mygene.MyGeneInfo()

	# Informations Recovery
	chro= line.split("\t")[1]
	start=int(line.split("\t")[2])
	end=int(line.split("\t")[3])
	region=line.split("\t")[7].split("(")[0]

	# If editing sites region is not Intergenic
	if region != "Intergenic" and region != "NA":

		# Edited Gene recovery
		tmp=("").join(line.split("\t")[7].split("(")[1])
	
		if "," in tmp:
			ncbi_gene= tmp.split(",")[0]
		elif ")" in tmp:
		
			ncbi_gene= tmp.split(")")[0]
	else:
		ncbi_gene="-"

	#If TSS associated gene is the same as annotated gene
	if ncbi_gene==line.split("\t")[10]:
		refgene=line.split("\t")[15]
		refgene_ok=True		
		
		
	# If editing sites is not annotated on Intergenic region
	
	if ncbi_gene !="-":
		if not refgene_ok:

			# NCBI to refgene conversion
			tab=[]
			refgene=""
			gene=mg.query(ncbi_gene)
			tab=(gene['hits'])

			if tab:
				dic=tab[0]
				refgene=dic['symbol']
			else:
			# if not refgene and NCBI equivalence Writing NCBI reference
				refgene=ncbi_gene
	else:
		refgene	="-"

	list_info.append(chro)
	list_info.append(start)
	list_info.append(end)
	list_info.append(refgene)
	list_info.append(region)
	return list_info
#########################################################################################################################
# Edited gene recovery and NCBI to refgene reference conversion

def conversion(f1,Name,path):

	path_out=path+"All/"
	if not os.path.isdir(path_out):
		os.mkdir(path_out)

	print ("Genes Conversion: this stage can take a long time")
	print("Please wait")

	mg = mygene.MyGeneInfo()
	f_out=open(path_out+Name+".bed","w")

	with open (f1,"r") as file1:
		#for each Editing Sites annotated
		for line in file1:
		
			if not "Gene" in line:
			
				list_info=ncbi_to_refgene(line)

				# Writing on file out: chromosome position position refgene region of each Editing Sites annotated
				f_out.write("%s\t%d\t%d\t%s\t%s\n"%(list_info[0],list_info[1],list_info[2],list_info[3],list_info[4]))

	f_out.close()

	# Output file sorting
	cmd = "cat "+path_out+Name+".bed|sort -k1,1 -k2,2n>"+path_out+Name+"_annotation.bed"
	cmd1 = "rm "+path_out+Name+".bed"
	subprocess.call(cmd,shell=True)
	subprocess.call(cmd1,shell=True)



###################################################################################################################
if __name__ == "__main__":
	
	f1="/media/bioinf/faef2f7b-527b-4b8a-859e-72ad73bd98bc/FINAL/Results/Known/Editing_sites/Annotation/Editing_sites_Na_Annotation.csv"
	name="Na"
	path_annotation="/media/bioinf/faef2f7b-527b-4b8a-859e-72ad73bd98bc/TEST/"
	conversion(f1,name,path_annotation)



