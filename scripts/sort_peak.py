#!/usr/bin/python 
# -*- coding: utf-8 -*-


###########################################################
#Mandy Cadix
#
#usage: sort_peak.py [-h] -input <INPUT_FILE> -output <OUTPUT_FILE> 
#
#File improvement (bed format)
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -input <INPUT_FILE>, --input <INPUT_FILE>
#                        Input file (bed format).
#  -output <OUTPUT_FILE>, --outputFile <OUTPUT_FILE>
#                        Output file (path/to/file/nameFile.bed).
#
###########################################################

import argparse
import re

##Tested with python 2.7.13
sort_peakVersion = '1.1 - 2016'

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='SCRIPT TO ORGANIZE COLUMNS OF GENE FILE', epilog='Version '+sort_peakVersion)
	parser.add_argument("-i","--input", type=str, action="store", default="", help="input file (bed file)",required=True)
	parser.add_argument("-o", "--output", type=str, action="store", default="", help ="output file - bed file with all information on peaks",required=True)

	dargs = vars(parser.parse_args())

	input_file=open(dargs["input_file"],"r")
	output_file=open(dargs["output_file"], "w")

	for line1 in input_file:
		tabline=line1.split("\t")
		nbCol=len(tabline)
		nbSamples=nbCol-9

	input_file.close()

	input_file=open(dargs["input_file"],"r")

	for line in input_file:
		tabline=line.split("\t")
		output_file.write(tabline[0]+"\t"+tabline[1]+"\t"+tabline[2]+"\t"+tabline[3]+"\t"+tabline[4]+"\t"+tabline[5]+"\t"+tabline[6]+"\t"+tabline[8]+"\t")
	
		for i in range(nbSamples):
			output_file.write(tabline[9+i].replace("\n","")+"\t")
		output_file.write(tabline[7]+"\n")
	
	input_file.close()
	output_file.close()
