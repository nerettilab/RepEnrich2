#!/usr/bin/env python
import argparse
import csv
import os
import shlex
import subprocess
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser(description='Part I: Prepartion of repetive element psuedogenomes and repetive element bamfiles.  This script prepares the annotation used by downstream applications to analyze for repetitive element enrichment. For this script to run properly bowtie must be loaded.  The repeat element psuedogenomes are prepared in order to analyze reads that map to multiple locations of the genome.  The repeat element bamfiles are prepared in order to use a region sorter to analyze reads that map to a single location of the genome.You will 1) annotation_file: The repetitive element annotation file downloaded from RepeatMasker.org database for your organism of interest. 2) genomefasta: Your genome of interest in fasta format, 3)setup_folder: a folder to contain repeat element setup files  command-line usage EXAMPLE: python master_setup.py /users/nneretti/data/annotation/mm9/mm9_repeatmasker.txt /users/nneretti/data/annotation/mm9/mm9.fa /users/nneretti/data/annotation/mm9/setup_folder', prog='getargs_genome_maker.py')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('annotation_file', action= 'store', metavar='annotation_file', help='List annotation file. The annotation file contains the repeat masker annotation for the genome of interest and may be downloaded at RepeatMasker.org  Example /data/annotation/mm9/mm9.fa.out')
parser.add_argument('genomefasta', action= 'store', metavar='genomefasta', help='File name and path for genome of interest in fasta format.  Example /data/annotation/mm9/mm9.fa')
parser.add_argument('setup_folder', action= 'store', metavar='setup_folder', help='List folder to contain bamfiles for repeats and repeat element psuedogenomes. Example /data/annotation/mm9/setup')
parser.add_argument('--threads', action='store', dest='threads', metavar='threads', default='1', type=int, help='Number of cores to be used for index building if possible.')
parser.add_argument('--nfragmentsfile1', action= 'store', dest='nfragmentsfile1', metavar='nfragmentsfile1', default='./repnames_nfragments.txt', help='Output location of a description file that saves the number of fragments processed per repname.  Default ./repnames_nfragments.txt')
parser.add_argument('--gaplength', action= 'store', dest='gaplength', metavar='gaplength', default= '200', type=int, help='Length of the spacer used to build repeat psuedogeneomes.  Default 200')
parser.add_argument('--flankinglength', action= 'store', dest='flankinglength', metavar='flankinglength', default= '25', type=int, help='Length of the flanking region adjacent to the repeat element that is used to build repeat psuedogeneomes.  The flanking length should be set according to the length of your reads.  Default 25')
parser.add_argument('--is_bed', action= 'store', dest='is_bed', metavar='is_bed', default= 'FALSE', help='Is the annotation file a bed file.  This is also a compatible format.  The file needs to be a tab seperated bed with optional fields.  Ex. format chr\tstart\tend\tName_element\tclass\tfamily.  The class and family should identical to name_element if not applicable.  Default FALSE change to TRUE')
args = parser.parse_args()

# parameters and paths specified in args_parse
gapl = args.gaplength
flankingl = args.flankinglength 
annotation_file = args.annotation_file
genomefasta = args.genomefasta
setup_folder = args.setup_folder
threads = args.threads
nfragmentsfile1 = args.nfragmentsfile1
is_bed = args.is_bed

################################################################################
# check that the programs we need are available
try:
    subprocess.call(shlex.split("bowtie2 --version"), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
except OSError:
    print "Error: Bowtie2 or BEDTools not loaded"
    raise

################################################################################
# Define a text importer
csv.field_size_limit(sys.maxsize)
def import_text(filename, separator):
    for line in csv.reader(open(os.path.realpath(filename)), delimiter=separator, 
                           skipinitialspace=True):
        if line:
            yield line
# Make a setup folder
if not os.path.exists(setup_folder):
	os.makedirs(setup_folder)

################################################################################
# load genome into dictionary
print "loading genome..."
g = SeqIO.to_dict(SeqIO.parse(genomefasta, "fasta"))

print "Precomputing length of all chromosomes..."
idxgenome = {}
lgenome = {}
genome = {}
allchrs = g.keys()
k = 0
for chr in allchrs:
    genome[chr] = str(g[chr].seq)
    del g[chr]
    lgenome[chr] = len(genome[chr])
    idxgenome[chr] = k
    k = k + 1
del g

################################################################################
# Build a bedfile of repeatcoordinates to use by RepEnrich region_sorter
if is_bed == "FALSE":
	repeat_elements= []
	fout = open(os.path.realpath(setup_folder + os.path.sep + 'repnames.bed'), 'w')
	fin = import_text(annotation_file, ' ')
	x = 0
	rep_chr = {}
	rep_start = {}
	rep_end = {}
	x = 0
	for line in fin:
	    if x>1:
	        line9 = line[9].replace("(","_").replace(")","_").replace("/","_")
	        repname = line9
		if not repname in repeat_elements:
			repeat_elements.append(repname)
	        repchr = line[4]
	        repstart = int(line[5])
	        repend = int(line[6])
		print >> fout, str(repchr) + '\t'+str(repstart)+ '\t'+str(repend)+ '\t'+str(repname)
	        if rep_chr.has_key(repname):
	            rep_chr[repname].append(repchr)
	            rep_start[repname].append(int(repstart))
	            rep_end[repname].append(int(repend))
	        else:
	            rep_chr[repname] = [repchr]
	            rep_start[repname] = [int(repstart)]
	            rep_end[repname] = [int(repend)]
	    x +=1
if is_bed == "TRUE":
	repeat_elements= []
	fout = open(os.path.realpath(setup_folder + os.path.sep + 'repnames.bed'), 'w')
	fin = open(os.path.realpath(annotation_file), 'r')
	x =0
	rep_chr = {}
	rep_start = {}
	rep_end = {}
	x =0
	for line in fin:
		line=line.strip('\n')
		line=line.split('\t')
	        line3 = line[3].replace("(","_").replace(")","_").replace("/","_")
	        repname = line3
		if not repname in repeat_elements:
			repeat_elements.append(repname)
	        repchr = line[0]
	        repstart = int(line[1])
	        repend = int(line[2])
		print >> fout, str(repchr) + '\t'+str(repstart)+ '\t'+str(repend)+ '\t'+str(repname)
	        if rep_chr.has_key(repname):
	            rep_chr[repname].append(repchr)
	            rep_start[repname].append(int(repstart))
	            rep_end[repname].append(int(repend))
	        else:
	            rep_chr[repname] = [repchr]
	            rep_start[repname] = [int(repstart)]
	            rep_end[repname] = [int(repend)]

fin.close()
fout.close()
repeat_elements = sorted(repeat_elements)
print "Writing a key for all repeats..."
#print to fout the binary key that contains each repeat type with the associated binary number; sort the binary key:
fout = open(os.path.realpath(setup_folder + os.path.sep + 'repgenomes_key.txt'), 'w')
x = 0
for repeat in repeat_elements:
    print >> fout, str(repeat) + '\t' + str(x)
    x +=1
fout.close()
################################################################################
# generate spacer for psuedogenomes
spacer = ""
for i in range(gapl):
    spacer = spacer + "N"

# save file with number of fragments processed per repname
print "Saving number of fragments processed per repname to " + nfragmentsfile1
fout1 = open(os.path.realpath(nfragmentsfile1),"w")
for repname in rep_chr.keys():
	rep_chr_current = rep_chr[repname]
	print >>fout1, str(len(rep_chr[repname])) + "\t" + repname
fout1.close()

# generate metagenomes and save them to FASTA files
k = 1
nrepgenomes = len(rep_chr.keys())
for repname in rep_chr.keys():
    metagenome = ""
    newname = repname.replace("(","_").replace(")","_").replace("/","_")
    print "processing repgenome " + newname + ".fa" + " (" + str(k) + " of " + str(nrepgenomes) + ")"
    rep_chr_current = rep_chr[repname]
    rep_start_current = rep_start[repname]
    rep_end_current = rep_end[repname]
    print "-------> " + str(len(rep_chr[repname])) + " fragments"
    for i in range(len(rep_chr[repname])):
        try:
            chr = rep_chr_current[i]
            rstart = max(rep_start_current[i] - flankingl, 0)
            rend = min(rep_end_current[i] + flankingl, lgenome[chr]-1)
            metagenome = metagenome + spacer + genome[chr][rstart:(rend+1)]
        except KeyError:
            print "Unrecognised Chromosome: "+chr
            pass
    
	# Convert metagenome to SeqRecord object (required by SeqIO.write)
    record = SeqRecord(Seq(metagenome, IUPAC.unambiguous_dna), id = "repname", name = "", description = "")
    print "saving repgenome " + newname + ".fa" + " (" + str(k) + " of " + str(nrepgenomes) + ")"
    fastafilename = os.path.realpath(setup_folder + os.path.sep + newname + ".fa")
    SeqIO.write(record, fastafilename, "fasta")
    print "indexing repgenome " + newname + ".fa" + " (" + str(k) + " of " + str(nrepgenomes) + ")"
    #perform bowtie2 indexing
    bt2filename = os.path.realpath(setup_folder + os.path.sep + newname)
    command = shlex.split('bowtie2-build -f --threads ' + str(threads) + ' ' + fastafilename + ' ' + bt2filename)
    #######
    p = subprocess.Popen(command).communicate()
    k += 1

print "... Done"

