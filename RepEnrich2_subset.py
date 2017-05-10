#!/usr/bin/env python
import argparse
import csv
import numpy
import os
import shlex
import shutil
import subprocess
import sys

parser = argparse.ArgumentParser(description='Prepartion for downstream running of RepEnrich2 by subsetting reads to uniquely mapped and multi-mapped. For this script to run properly samtools and bedtools must be loaded.  You will need 1) An alignment file. 2) A defined MAPQ threshold (30 by default)   command-line usage EXAMPLE: python RepEnrich_subset.py /example/path/data/alignment.sam 30', prog='RepEnrich2_subset.py')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('alignment_file', action= 'store', metavar='alignment_file', help='bam file of aligned reads. Example /example/path/alignment.bam')
parser.add_argument('MAPQ', action= 'store', metavar='MAPQ', default='30', help='MAPQ threshold for subsetting uniquely mapping and multi-mapping reads.  This will be aligner-dependent as many aligners treat the MAPQ field differently (Default = 30) ')
parser.add_argument('sample_name', action='store', metavar='sample_name', help='The name of your sample (used to generate the name of the output files)')
parser.add_argument('--pairedend', action= 'store', dest='pairedend', default= 'FALSE', help='Designate this option for paired-end data. Default FALSE change to TRUE')

args = parser.parse_args()


# parameters specified in args_parse
alignment_file = args.alignment_file
MAPQ = args.MAPQ
sample_name = args.sample_name
pairedend = args.pairedend

################################################################################
# check for loaded dependencies
try:
	subprocess.call(shlex.split("samtools --version"), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
except OSError:
	print "Error: Samtools not loaded"
	raise
try:
	subprocess.call(shlex.split("bedtools --version"), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
except OSError:
	print "Error: Bedtools not loaded"
	raise
################################################################################


#output .bam (unique) and .fastq/s (multimapping) alignments

if pairedend == "FALSE":
	print "Subsetting uniquely mapping reads above MAPQ threshold..."
	fileout= sample_name + '_unique.bam'
	with open(fileout, 'w') as stdout:
		#filter for mapped reads, with MAPQ above the specified value
		command = shlex.split("samtools view -F 4 -bq " + MAPQ + " " + alignment_file)
		p = subprocess.Popen(command, stdout=stdout)
	p.communicate()
	stdout.close()
	print "Done."


if pairedend == "TRUE"	:
	print "Subsetting uniquely mapping paired-end reads in proper pair above MAPQ threshold..."
	fileout= sample_name + '_unique.bam'
	with open(fileout, 'w') as stdout:
		#filter for paired reads with read mapped in proper pair, with MAPQ above the specified value
		command = shlex.split("samtools view -f 3 -bq " + MAPQ + " " + alignment_file)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()
	print "Done."


##############################
if pairedend == "FALSE":
	print "Subsetting multi-mapping reads..."
	fileout2= sample_name + '_map.bam'
	fileout3= sample_name + '_multimap_filtered.bam'
	fileout4= sample_name + '_multimap.fastq'
	with open(fileout2, 'w') as stdout:
		#filter for mapped reads only
		command = shlex.split("samtools view -h -F 4 -b " + alignment_file)
		p = subprocess.Popen(command, stdout=stdout)

		p.communicate()
		stdout.close()

	#filter from the mapped reads for everything that does not have a MAPQ meeting the threshold; everything should already map from the last filter
	#-U option requires samtools 1.3.0+
	command = shlex.split("samtools view -U " + fileout3 + " -bq " + MAPQ + " " + fileout2)
	#output to fileout3 and silence stdout
	p = subprocess.Popen(command, stdout=open(os.devnull, 'wb'))
	p.communicate()

	#convert bamfile to fastq for downstream use
	command = shlex.split("bedtools bamtofastq -i " + fileout3 + " -fq " + fileout4)
	p=subprocess.Popen(command)
	p.communicate()
	print "Done."

	print "Performing cleanup..."
	#cleanup intermediate files
	command = shlex.split("rm " + fileout2)
	p=subprocess.Popen(command)
	p.communicate()
	stdout.close()
	command = shlex.split("rm " + fileout3)
	p.communicate()
	stdout.close()
	print "Done."

##############################
if pairedend == "TRUE":
	print "Subsetting paired-end multi-mapping reads..."
	fileout2= sample_name + '_map.bam'
	fileout3= sample_name + '_multimap_filtered.bam'
	fileout4= sample_name + '_multimap_sorted.bam'
	fileout5= sample_name + '_multimap_R1.fastq'
	fileout6= sample_name + '_multimap_R2.fastq'

	with open(fileout2, 'w') as stdout:
    #filter for paired / mapped in proper pair only
		command = shlex.split("samtools view -f 3 -b " + alignment_file)
		p = subprocess.Popen(command, stdout=stdout)
		p.communicate()
		stdout.close()

	#filter from the mapped reads for everything that does not have a MAPQ meeting the threshold for being uniquely mapped
	command = shlex.split("samtools view -U " + fileout3 + " -bq " + MAPQ + " " + fileout2)
	#output to fileout3 and silence stdout
	p = subprocess.Popen(command, stdout=open(os.devnull, 'wb'))
	p.communicate()
	stdout.close()

	print "Performing sort on bamfile..."
	#sort bamfile by name for use with bamtofastq (sorting only needed for paired end reads)
	command = shlex.split("samtools sort -n " + fileout3 + " -o " + fileout4)
	p = subprocess.Popen(command)
	p.communicate()
	

	print "Converting to fastq..."
	#convert bamfile to fastq for downstream use
	command = shlex.split("bedtools bamtofastq -i " + fileout4 + " -fq " + fileout5 + " -fq2 " + fileout6)
	p=subprocess.Popen(command)
	p.communicate()
	

	print "Performing cleanup..."
	#cleanup intermediate files
	command = shlex.split("rm " + fileout2)
	p=subprocess.Popen(command)
	p.communicate()
	command = shlex.split("rm " + fileout3)
	p=subprocess.Popen(command)
	p.communicate()
	command = shlex.split("rm " + fileout4)
	p=subprocess.Popen(command)
	p.communicate()

	print "Done."


