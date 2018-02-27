#!/usr/bin/python 
# -*- coding: utf-8 -*-

###########################################################
#Mandy Cadix
#
#usage: sort_gene.py [-h] -i <INPUT_FILE> -o <OUTPUT_FILE>
#
#Add & sort the gene file
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -i <INPUT_FILE>, --inputFile <INPUT_FILE>
#                        Input file (bed format).
#  -o <OUTPUT_FILE>, --outputFile <OUTPUT_FILE>
#                        Output file (path/to/file/nameFile.bed).
#
###########################################################

import argparse
import os
import re

## Tested with python 2.7.12
sort_geneVersion = '1.1 - 01/10/2015'

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'SCRIPT TO ORGANIZE COLUMNS OF GENE FILE', epilog='Version '+sort_geneVersion)

	parser.add_argument("-i","--input", type = str, action = "store", default = "", help = "input file (bed file)", required = True)
	parser.add_argument("-o", "--output", type = str, action = "store", default = "", help = "output file - bed file with all information on genes", required = True)
	dargs = vars(parser.parse_args())

	input_file = open(dargs["input_file"],"r")
	output_file = open(dargs["output_file"], "w")

	for line1 in input_file:
		tabline = line1.split("\t")
		nbCol = len(tabline)
		nbSamples = nbCol-7

	input_file.close()

	input_file = open(dargs["input_file"],"r")

	for line in input_file:
		tabline = line.split("\t")
		output_file.write(tabline[1] +"\t"+ tabline[2] +"\t"+ tabline[3] +"\t"+ tabline[5] +"\t"+ tabline[5] +"\t"+ tabline[4] +"\t"+ tabline[0] +"\t"+ tabline[nbCol-1].replace("\n","\t"))
		for i in range(nbSamples):
			output_file.write(tabline[6+i] +"\t")
		tabName = line.split("_")
		output_file.write(tabName[0] +"_"+ tabName[1] +"_"+ tabName[2] +"_"+ tabName[3] +"\n")

	input_file.close()
	output_file.close()
