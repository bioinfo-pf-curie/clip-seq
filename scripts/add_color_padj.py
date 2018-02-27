#! /usr/bin/python
# -*- coding: utf8 -*-

###########################################################
#
#usage: add_color_padj.py [-h] -input <INPUT_FILE> -output <OUTPUT_FILE>
#
#Add color code in function of gene (Down, Up, Nothing expressed), padj<0.05.
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -input <INPUT_FILE>, --inputFile <INPUT_FILE>
#                        Input file (chr start end peakName ...).
#  -output <OUTPUT_FILE>, --outputFile <OUTPUT_FILE>
#                        Output file (path/to/file/nameFile).
#
#
###########################################################

###########################################################
## Import
import argparse
import os

## Tested with python 2.7.5
addColorVersion = '1.0 - 14/12/2015'

if __name__ == '__main__':

	########### arguments ####################
	parser = argparse.ArgumentParser(description='Add color code', epilog='Version '+addColorVersion)
	
	parser.add_argument('-input', '--inputFile', type=str, action="store", default="", help='Input file (chr start end peakName ...).', required=True)
	parser.add_argument('-output', '--outputFile', type=str, action="store", default="", help='Output file (path/to/file/nameFile).', required=True)

	dargs = vars(parser.parse_args())

	inputFile = open(dargs["inputFile"],"r")

	lineDescription = inputFile.readline()

	outputFile = open(dargs["outputFile"],"w")
	outputFile.write("peak"+"\t"+lineDescription.replace("\n","\t")+"rgbColor"+"\n")

	for line in inputFile :

		tabLine = line.split("\t")
		if (float(tabLine[10])<0.05 and float(tabLine[21])<0):
			outputFile.write(line.replace("\n","\t")+"0,0,255"+"\n")
		elif (float(tabLine[10])<0.05 and float(tabLine[21])>0):
			outputFile.write(line.replace("\n","\t")+"255,64,64"+"\n")
		else:
			outputFile.write(line.replace("\n","\t")+"0,0,0"+"\n")

	outputFile.close()
	inputFile.close()


