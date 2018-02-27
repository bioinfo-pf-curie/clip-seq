#! /usr/bin/python
# -*- coding: utf8 -*-

###########################################################
#Mandy Cadix
#
#usage: extended_gene.py [-h] -inputA <INPUT_FILE_A> -inputB <INPUT_FILE_B> -output <OUTPUT_FILE>
#
#Treat genes which have an extension many pb.
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -inputA <INPUT_FILE_A>, --inputFileA <INPUT_FILE_A>
#                        Input first file (chr start end geneName nbFusion strand).
#  -inputB <INPUT_FILE_B>, --inputFileB <INPUT_FILE_B>
#                        Input second file (chr start end geneName nbFusion strand).
#  -output <OUTPUT_FILE>, --outputFile <OUTPUT_FILE>
#                        Output file (path/to/file/nameFile).
#
#
###########################################################

import argparse
import os

## Tested with python 2.7.5
extended_geneVersion = '1.0 - 26/10/2015'

if __name__ == '__main__':

	########### arguments ####################
	parser = argparse.ArgumentParser(description = 'Treat genes which have an extension many pb', epilog = 'Version ' + extended_geneVersion)
	
	parser.add_argument('-inputA', '--inputFileA', type = str, action = "store", default = "", help = 'Input first file (whole gene)', required = True)
	parser.add_argument('-inputB', '--inputFileB', type = str, action = "store", default = "", help = 'Input second file (whole gene with new coordinate)', required = True)
	parser.add_argument('-output', '--outputFile', type = str, action = "store", default = "", help = 'Output file (path/to/file/nameFile)', required = True)

	dargs = vars(parser.parse_args())

	inputFileA = open(dargs["inputFileA"], "r")
	inputFileB = open(dargs["inputFileB"], "r")
	outputFile = open(dargs["outputFile"], "w")

	for lineA in inputFileA :
		tabLineA = lineA.split("\t")
		inputFileB = open(dargs["inputFileB"], "r")
		for lineB in inputFileB :
			tabLineB = lineB.split("\t")
			if (tabLineA[3] == tabLineB[3]):
				if ("+" in tabLineA[3]):
					outputFile.write(tabLineA[0] + "\t" + tabLineA[1] + "\t" + tabLineB[2] + "\t" + tabLineA[3] + "\t" + tabLineA[4] + "\t" + tabLineA[5])
				else:
					outputFile.write(tabLineA[0] + "\t" + tabLineB[1] + "\t" + tabLineA[2] + "\t" + tabLineA[3] + "\t" + tabLineA[4] + "\t" + tabLineA[5])



	inputFileA.close()
	inputFileB.close()
	outputFile.close()

