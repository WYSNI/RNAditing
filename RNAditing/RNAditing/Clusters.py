#!/bin/python
#coding:utf_8

from collections import defaultdict
import numpy as np
from random import shuffle
import os,subprocess,mygene,glob
import Genes_recovery as Genes

def Clusters_Annotation(file_tmp,popu,path_cluster):
	
	print("*** Annotation ***\n")

	outFile=path_cluster+"Tmp_Editing_islands_"+popu+".bed"
	file_ann=path_cluster+"Editing_islands_"+popu+"_Annotation.csv"
	Tmp_ann=path_cluster+"Tmp_ann.csv"
	
	cmd ="annotatePeaks.pl "+file_tmp+ " hg19 -annStats "+path_cluster+"Annotation_"+popu+"_Stats.csv >"+Tmp_ann+" 2>"+path_cluster+"trash.txt"
	subprocess.call(cmd,shell=True)

	cmd="rm "+path_cluster+"trash.txt"
	subprocess.call(cmd, shell=True)

	print("\t- Annotation for %s file OK: %s and %sAnnotation_%s_Stats.csv created\n\n"%(popu,file_ann,path_cluster,popu))
	
	cmd ="cat "+Tmp_ann+"|sort -k1,1n>"+file_ann
	subprocess.call(cmd,shell=True)
	
	cmd ="rm "+Tmp_ann
	subprocess.call(cmd,shell=True)
	
	# Edited genes recovery

	print("*** Genes Conversion: this stage can take a long time ***")
	print("Please wait\n")

	f_out=open(outFile,"w")

	with open (file_ann,"r") as file1:

		#for each Editing Sites annotated
		for line in file1:
			if not "Gene" in line:
			
				# Informations Recovery
				island=line.split("\t")[0]
				list_info=Genes.ncbi_to_refgene(line)

				# Writing on file out: chromosome position position refgene region of each Editing Sites annotated
				f_out.write("%s\t%d\t%d\t%s\t%s\t%s\t\n"%(list_info[0],list_info[1],list_info[2],island,list_info[3],list_info[4]))
	
	return outFile
#################################################################################################

#Dictionnary creation {"chr1":[pos,pos,pos,pos],"chr12":[pos,pos,pos,pos,pos,pos]...}
def create_Dic_Pos(file_in):
	dic_pos={}
	chro=""
	List_pos=[]
	i=0
	nb_line=0
	with open (file_in,"r") as bedfile:
		#Number of lines in file
		for l in bedfile:
			nb_line+=1
		bedfile.seek(0,0)

		#For each line
		for line in bedfile:
			i+=1
			l_splited=line.split()
			
			if chro == "":
				chro=l_splited[0]
			elif chro != l_splited[0]:
				dic_pos[chro]=List_pos
				chro=l_splited[0]
				List_pos=[]
			List_pos.append(int(l_splited[1]))
			
			#if current line is the last one of the file
			if i==nb_line:
				dic_pos[chro]=List_pos
	return (dic_pos) 

#################################################################################################

def printClusters(clusterDict,path_out,popu):
	path_clusters=path_out+"Clusters/"
	#path_clusters=path_out+"Clusters_pvalue/"
	if not os.path.isdir(path_clusters):
		os.mkdir(path_clusters)

	outFile=path_clusters+"Editing_islands_"+popu+".bed"
	outTmp=path_clusters+"Tmp.bed"
	outTmp_sorted=path_clusters+"Tmp_sorted.bed"
	
	out_Tmp=open(outTmp,"w")
	
	for chr,Dic_clusters in zip(clusterDict.keys(),clusterDict.values()):
	
		for key,clusters in zip(Dic_clusters.keys(),Dic_clusters.values()):

			end=max(clusters)
			start=min(clusters)
			
			length = end - start
			editingRate=float(len(clusters))/float(length)
			
			out_Tmp.write("\t".join([str(chr),str(start),str(end),str(key),str(length),str(len(clusters)),'%1.2f'%float(editingRate),"\n"]))

	out_Tmp.close()
	
	cmd="cat "+outTmp+"|sort -k4,4n>"+outTmp_sorted
	subprocess.call(cmd,shell=True)

	#Annotation

	Tmp_clusters=Clusters_Annotation(outTmp_sorted,popu,path_clusters)

	#Add further informations

	Clusters_file=path_clusters+"Editing_islands_"+popu+".bed"
	Clusters=open(Clusters_file,"w")
	Clusters.write("\t".join(["#Chr","Start","Stop","IslandID","Gene Symbol","Region","Cluster Length","Number of Editing_sites","Editing_rate","\n"]))
	
	with open (Tmp_clusters,"r") as f1:
		
		with open (outTmp_sorted,"r") as f2:
			lines1=f1.readlines()
			lines2=f2.readlines()
			
			for linef1 in lines1:
				for linef2 in lines2:
					if str(linef1.split("\t")[3]) == str(linef2.split("\t")[3]):
						
						Clusters.write(("\t").join(linef2.split("\t")[0:3])+"\t"+("\t").join(linef1.split("\t")[3:6])+"\t"+("\t").join(linef2.split("\t")[4:]))
						break
				
	Clusters.close()

	cmd ="rm "+Tmp_clusters+" "+outTmp_sorted+" "+outTmp
	subprocess.call(cmd, shell= True)
	
	return Clusters_file

#################################################################################################

def createClusters(Dic_pos,eps=50,minSamples=5):
	islandCounter=0
	Dic_clustersByChr={}

	# For each Chromosome: 
	for chr in Dic_pos.keys():
		clusterDict={}
		posList = Dic_pos[chr] #position of all variants from that chromosome
		labels = getLabels(posList,eps,minSamples) #actually doing db clustering
		n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
		
		if n_clusters_ > 0:
			#loop over labels and variants
			tmpDict=defaultdict(list)
			for var,label in zip(Dic_pos[chr],labels):
				if label==-1:
					continue
		            
				tmpDict[label].append(var)
		            
			#set new label for clusterdict, to avoid overwriting
			for label in tmpDict.keys():
				clusterDict[islandCounter]=tmpDict.pop(label)
				islandCounter+=1
			Dic_clustersByChr[chr]=clusterDict
	
	return Dic_clustersByChr
	
#################################################################################################

def getLabels(positionList,eps=10, minSamples=5):
	"""Perform DBSCAN clustering from vector array.

	Parameters
	----------
	X: array [int1,int1]
	    Array of Samples. 
	    In this case it should be the positions of the variations in the genome per chromosome

	eps: float, optional
	    The maximum distance between two samples for them to be considered
	    as in the same neighborhood.

	minSamples: int, optional
	    The number of samples in a neighborhood for a point to be considered
	    as a core point.


	Returns
	-------
	core_samples: array [n_core_samples]
	    Indices of core samples.

	labels : array [n_samples]
	    Cluster labels for each point.  Noisy samples are given the label -1.

	"""
	if not eps > 0.0:
	    raise ValueError("eps must be positive.")

	X = np.asarray(positionList)	
	n = X.shape[0] #get number of elements
	index_order=range(n)# Index creation

	shuffle(index_order)#random shuffle index order


	distanceMatrix = calculate1dDistanceMatrix(X,eps)
	
	# Calculate neighborhood for all samples. This leaves the original point
	# in, which needs to be considered later (i.e. point i is the
	# neighborhood of point i. While True, its useless information)

	#distanceMatrix = [np.where(x <= eps)[0] for x in distanceMatrix]

	# Initially, all samples are noise.(-1)
	labels = -np.ones(n, dtype=np.int)

	# A list of all core samples found.
	core_samples = []

	# label_num is the label given to the new cluster
	label_num = 0

	# Look at all samples and determine if they are core.
	# If they are then build a new cluster from them.
	nb=0

	for index in index_order:
		
		# Already classified
		if labels[index] != -1: 
			continue

		# get neighbors from distanceMatrix or ballTree
		index_neighborhood = []

		index_neighborhood = distanceMatrix[index]# Position neighborhood recovery 


		# Too few samples to be core
		if len(index_neighborhood) < minSamples:
			continue

		core_samples.append(index)
		
		labels[index] = label_num
		
		# candidates for new core samples in the cluster.
		candidates = [index]
		
		while len(candidates) > 0:
			
			new_candidates = []
			# A candidate is a core point in the current cluster that has
			# not yet been used to expand the current cluster.
			for c in candidates:
				
				c_neighborhood = []

				c_neighborhood = distanceMatrix[c]
				noise = np.where(labels[c_neighborhood] == -1)[0] #indexes of candidate neigbours which do not belong to a cluster yet
				noise = c_neighborhood[noise]
				labels[noise] = label_num
				for neighbor in noise:
					n_neighborhood = []

					n_neighborhood = distanceMatrix[neighbor]

					# check if its a core point as well
					if len(n_neighborhood) >= minSamples:

						# is new core point
						new_candidates.append(neighbor)
						core_samples.append(neighbor)

			# Update candidates for next round of cluster expansion.
			candidates = new_candidates

		# Current cluster finished.
		# Next core point found will start a new cluster.
		label_num += 1
	#return core_samples, labels
        
	return labels

################################################################################################# 
   
def calculate1dDistanceMatrix(lst,eps):
	'''
	creates a distance matrix for the given vector
	:param lst: vector of samples       
	:return: np.array(diffMatrix)
	'''
	if not isinstance(lst, (list, tuple, np.ndarray)): #Is it a List or a Tuple?
		raise TypeError("Paramer has to be eithe a List or a Tuple found %s" % type(lst))
	if not all(isinstance(item, (int,float)) for item in lst): # Is it a list of number?
		raise TypeError("List should only contain numbers")
	lst = np.asarray(lst)
	diffMatrix=[]

	#For each position
	for l1 in lst:

		diffList=[]
		diffList= abs(lst-l1) # Check the distance beetween each position of the list an the current position
		diffList = np.where(diffList<=eps)[0] # Each distance lower than eps are recovered on an array
		diffMatrix.append(diffList) # The diffList array is recovered on an array
	return np.asarray(diffMatrix)

#################################################################################################
def Cluster_recovery(path_results):
	#list_files=glob.glob(path_results+"Editing_sites/Editing*.bed")
	list_files=glob.glob(path_results+"Editing_sites/Editing*.bed")
	for file in list_files:
		if "RNAEditor_" in file:
			popu=file[:-4].split("RNAEditor_")[1]

		elif "De_novo_" in file:
			popu=file[:-4].split("De_novo_")[1]
			
		else:
			popu=file[:-4].split("Editing_Sites_")[1]
			#popu=file[:-4].split("_")[-1]

		print("\n##### Clusters Recovery for %s population ##### \n"%popu)
		
		
		Dic_pos=create_Dic_Pos(file)
		clusterDict=createClusters(Dic_pos,eps=50,minSamples=5)
		Clusters_file=printClusters(clusterDict,path_results,popu)
		
		print("\t- Clusters recovery for %s population OK: %s created\n"%(popu,Clusters_file))
		
#################################################################################################

if __name__=="__main__":
	#file_in="/media/bioinf/faef2f7b-527b-4b8a-859e-72ad73bd98bc/TEST/Results/Editing_sites/Editing_Sites_REDItools_Na_RNAEditor_Na.bed"
	path="/media/bioinf/faef2f7b-527b-4b8a-859e-72ad73bd98bc/TEST/RNAEditor_REDItools/"
	Cluster_recovery(path)
