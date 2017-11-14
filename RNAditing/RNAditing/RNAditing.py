#!/bin/python2.7
# coding:utf_8

import sys
import RNAEditor 
import Editing_sites_knowns as ESK
import All_analysis as All
import usage

# VERSION: 1.0.0
# Author: Benjamin Delaune Bioinformatic student team 11 CRCINA
# Date : 30/08/17

########################################################################################################
# Pipeline feature

def ref():
	print """
	########################################################################
	# RNAditing: RNA Editing sites research and prediction on RNA-Seq data #
	# VERSION: 1.0.0						       #
	# Author: Benjamin Delaune, Bioinformatic Student, Team 11 CRCINA      #
	######################################################################## 
	
	"""

########################################################################################################
# MAIN 

def main():
	# Pipeline feature
	ref()

	#Argument number verification
	if len(sys.argv) <2:
		usage.RNAditing_usage()
		sys.exit(2)
	
	#Option choice
	if sys.argv[1] == "-R" or sys.argv[1] == "--RnaEditor":
		RNAEditor.startAnalysis(sys.argv[2:])

	elif sys.argv[1] == "-E" or sys.argv[1] == "--ESknowns":
		ESK.startAnalysis(sys.argv[2:])

	elif sys.argv[1] == "-A" or sys.argv[1] == "--All":
		All.startAnalysis(sys.argv[2:])
	else:
		usage.RNAditing_usage()
		sys.exit(2)

	print("** Program complete **")

#########################################################################################################

if __name__=="__main__":
	main()













