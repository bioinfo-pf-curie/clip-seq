#!/usr/bin/python
# -*- coding: utf-8 -*-

###########################################################
#Mandy Cadix
#
#usage: selectFile.py [-h] -input <INPUT_FILE> -output <OUTPUT_FILE> 
#
#Select sequences (fasta.gz format) without problem problem of sequencing.
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -input <INPUT_FILE>, --input <INPUT_FILE>
#                        Input file (fastq.gz format).
#  -output <OUTPUT_FILE>, --outputFile <OUTPUT_FILE>
#                        Output file (path/to/file/nameFile.fastq).
#
###########################################################


import argparse
import gzip
import re

##Tested with python 2.7.9
selectFileVersion = '1.1 - 03/01/2018'

if __name__ == '__main__':
	##### arguments #####
	parser = argparse.ArgumentParser(description = 'SCRIPT TO SELECT THE SEQUENCES WITHOUT PROBLEM OF SEQUENCING', epilog='Version '+selectFileVersion)
	parser.add_argument('-input','--inputFile', type=str, action="store", default="", help = 'input file (fastq.gz format)', required = True)
	parser.add_argument('-output', '--outputFile', type=str, action="store", default="", help = "fastq file", required = True) 
	
	dargs = vars(parser.parse_args())


	inputFile = gzip.open(dargs["inputFile"],"rb")
	outputFile = open(dargs["outputFile"],"w")

#	input_file=gzip.open(input_file,'rb')
#	output_file = open(output_file, 'w')


flag = 0
y_counter = 0
n_counter = 0

for line in inputFile:	

	if (" " in line) and (":N:" in line):
		flag = 1
		n_counter += 1

	if (" " in line) and (":Y:" in line):
		flag = 0
		y_counter += 1

	if flag == 1:
		outputFile.write(line)

print "Total reads: ", n_counter + y_counter
print "Discarded reads: ", y_counter

inputFile.close()
outputFile.close()
