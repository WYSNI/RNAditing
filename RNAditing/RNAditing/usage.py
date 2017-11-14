#!bin/python
#coding:utf_8


# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

#####################################################################################################################################
# Help for use the full pipeline

def All_Analysis_usage():
	print """
	USAGE: RNAditing --All -a [options]
	Options:

	-a or --analysis		Use this option for use RNAEditor in commande line (needed)

		-c or --configure		Use this option with the -n option for give to RNAEditor the configuration file needed (needed)
		-f or --fastq			Used for specify the fastq folder used (needed)
		-s or --sample			Used for give the sample file path (needed)
		-o or --output			Used for give the future Bam files folder path used by REDItools (needed)
		-p or --percentage		Minimum percentage of a site for it to be considered representative of a population : 0-1 (needed)
		-e or --end			Paired-end or single-end sequencing. Write paired or single (needed)
	"""

######################################################################################################################################
# Help for use REDItools software

def REDItools_usage():
	print """
	USAGE: RNAditing --ESknowns [-m or -a or -n] [options]
	Options:

	-a or --all			Use this option for use REDItools in commande line and make a comparison by population

		-d or --Database		Database with editing sites already known (needed)
		-i or --Bam			Path to directory containing BAMs to analyzed (needed)
		-r or --Reference		Path to reference fasta file (needed)
		-s or --Splice			Path to splice file (splicesites.ss)
		-v or --variant			Path to variant file (dbSNP)
		-o or --output			Path to results folder (needed)
		-p or --percentage		Minimum percentage of a site for it to be considered representative of a population : 0-1 (needed)
		-c or --software_comp		Path to RNAEditor results (use if you want to make a softwares results comparison)
		
	-n or --analyse			Use this option for use REDItools in commande line

		-d or --Database		Database with editing sites already known (needed)
		-i or --Bam			Path to directory containing BAMs to analyzed (needed)
		-r or --Reference		Path to reference fasta file (needed)
		-s or --Splice			Path to splice file (splicesites.ss)
		-v or --variant			Path to variant file (dbSNP)
		-o or --output			Path to results folder (needed)

	-m or --merge			Use for make a comparison by population (For MN_3, MN_8, MN_9 for example)

		-d or --Database		Database with editing sites already known (needed)
		-o or --output			Path to results folder (needed)
		-p or --percentage		Minimum percentage of a site for it to be considered representative of a population : 0-1 (needed)

		### For the comparison with RNAEditor add :
		
		-c or --software_comp		Path to RNAEditor results (needed)
			
	"""
######################################################################################################################################
# Help for choose the options

def RNAditing_usage():
	print """
	USAGE: RNAditing -R or -E or -A [options]
	Options:

	-A or --All			To use all of the pipeline
	-R or --rnaEditor		Used for the prediction of RNA Editing sites
	-E or --ESknowns		Used for the discovery of RNA Editing sites already known on databases


	"""

######################################################################################################################################
# Help for use RNA Editor software

def RNAEditor_usage():
	print """
	USAGE: RNAditing --RnaEditor [-i or -m or -a or -n] [options]
	Options:

	-i or --GUI			Use this option alone for use the RNAEditor software with its interface.

	-a or --all			Use this option for use RNAEditor in commande line and make a comparison by population 

		-c or --configure		Use this option to give to RNAEditor the configuration file path (needed)
		-f or --fastq			Use for specify the FASTQ folder used (needed)	
		-p or --percentage		Minimum percentage of a site for it to be considered representative of a population : 0-1 (needed)
		-s or --sample			Used for give the samples file path (needed)
		-e or --end			Paired-end or single-end sequencing. Write paired or single (needed)

	-n or --analyse			Use this option for use RNAEditor in commande line  (needed)

		-c or --configure		Use this option to give to RNAEditor the configuration file path (needed)
		-f or --fastq			Use for specify the FASTQ folder used (needed)	
		-s or --sample			Used for give the samples file path (needed)
		-e or --end			Paired-end or single-end sequencing. Write paired or single (needed)


	-m or --merge			Use for make a comparison by population (MN_3, MN_8, MN_9 for example)

		-f or --fastq			Use for specify the FASTQ folder used (needed)	
		-p or --percentage		Minimum percentage of a site for it to be considered representative of a population : 0-1 (needed)
	


	"""
