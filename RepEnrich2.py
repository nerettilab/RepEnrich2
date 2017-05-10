#!/usr/bin/env python
import argparse
import csv
import numpy
import os
import shlex
import shutil
import subprocess
import sys

parser = argparse.ArgumentParser(description='Part II: Conducting the alignments to the psuedogenomes.  Before doing this step you will require 1) a bamfile of the unique alignments with index 2) a fastq file of the reads mapping to more than one location.  These files can be obtained using the following bowtie options [EXAMPLE: bowtie -S -m 1 --max multimap.fastq mm9 mate1_reads.fastq]  Once you have the unique alignment bamfile and the reads mapping to more than one location in a fastq file you can run this step.  EXAMPLE: python master_output.py /users/nneretti/data/annotation/hg19/hg19_repeatmasker.txt /users/nneretti/datasets/repeatmapping/POL3/Pol3_human/HeLa_InputChIPseq_Rep1 HeLa_InputChIPseq_Rep1 /users/nneretti/data/annotation/hg19/setup_folder HeLa_InputChIPseq_Rep1_multimap.fastq HeLa_InputChIPseq_Rep1.bam')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('annotation_file', action= 'store', metavar='annotation_file', help='List RepeatMasker.org annotation file for your organism.  The file may be downloaded from the RepeatMasker.org website.  Example: /data/annotation/hg19/hg19_repeatmasker.txt')
parser.add_argument('outputfolder', action= 'store', metavar='outputfolder', help='List folder to contain results. Example: /outputfolder')
parser.add_argument('outputprefix', action= 'store', metavar='outputprefix', help='Enter prefix name for data. Example: HeLa_InputChIPseq_Rep1')
parser.add_argument('setup_folder', action= 'store', metavar='setup_folder', help='List folder that contains the repeat element psuedogenomes. Example /data/annotation/hg19/setup_folder')
parser.add_argument('fastqfile', action= 'store', metavar='fastqfile', help='Enter file for the fastq reads that map to multiple locations. Example /data/multimap.fastq')
parser.add_argument('alignment_bam', action= 'store', metavar='alignment_bam', help='Enter bamfile output for reads that map uniquely. Example /bamfiles/old.bam')
parser.add_argument('--pairedend', action= 'store', dest='pairedend', default= 'FALSE', help='Designate this option for paired-end sequencing. Default FALSE change to TRUE')
parser.add_argument('--collapserepeat', action= 'store', dest='collapserepeat', metavar='collapserepeat', default= 'Simple_repeat', help='Designate this option to generate a collapsed repeat type.  Uncollapsed output is generated in addition to collapsed repeat type.  Simple_repeat is default to simplify downstream analysis.  You can change the default to another repeat name to collapse a seperate specific repeat instead or if the name of Simple_repeat is different for your organism.  Default Simple_repeat')
parser.add_argument('--fastqfile2', action= 'store', dest='fastqfile2', metavar='fastqfile2', default= 'none', help='Enter fastqfile2 when using paired-end option.  Default none')
parser.add_argument('--cpus', action= 'store', dest='cpus', metavar='cpus', default= "1", type=int, help='Enter available cpus per node.  The more cpus the faster RepEnrich performs. RepEnrich is designed to only work on one node. Default: "1"')
parser.add_argument('--allcountmethod', action= 'store', dest='allcountmethod', metavar='allcountmethod', default= "FALSE", help='By default the pipeline only outputs the fraction count method.  Consdidered to be the best way to count multimapped reads.  Changing this option will include the unique count method, a conservative count, and the total count method, a liberal counting strategy. Our evaluation of simulated data indicated fraction counting is best. Default = FALSE, change to TRUE')
parser.add_argument('--is_bed', action= 'store', dest='is_bed', metavar='is_bed', default= 'FALSE', help='Is the annotation file a bed file.  This is also a compatible format.  The file needs to be a tab seperated bed with optional fields.  Ex. format chr\tstart\tend\tName_element\tclass\tfamily.  The class and family should identical to name_element if not applicable.  Default FALSE change to TRUE')
args = parser.parse_args()

# parameters
annotation_file = args.annotation_file
outputfolder = args.outputfolder
outputfile_prefix = args.outputprefix
setup_folder = args.setup_folder
repeat_bed = setup_folder + os.path.sep + 'repnames.bed' 
unique_mapper_bam = args.alignment_bam
fastqfile_1 = args.fastqfile
fastqfile_2 = args.fastqfile2
cpus = args.cpus
b_opt = "-k 1 -p " +str(1) +" --quiet"
simple_repeat = args.collapserepeat
paired_end = args.pairedend
allcountmethod = args.allcountmethod
is_bed = args.is_bed

################################################################################
# check that the programs we need are available
try:
    subprocess.call(shlex.split("coverageBed -h"), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    subprocess.call(shlex.split("bowtie2 --version"), stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
except OSError:
    print "Error: Bowtie2 or BEDTools not loaded"
    raise

# check version of bedtools and store numbers in variables for determining syntax of bedCoverage
btv1 = subprocess.Popen('bedtools --version | cut -d" " -f2 | cut -d"v" -f2 | cut -d"." -f1', shell=True, stdout=subprocess.PIPE)
v1 = btv1.stdout.read()
btv2 = subprocess.Popen('bedtools --version | cut -d" " -f2 | cut -d"v" -f2 | cut -d"." -f2', shell=True, stdout=subprocess.PIPE)
v2 = btv2.stdout.read()

################################################################################
# define a csv reader that reads space deliminated files
print 'Preparing for analysis using RepEnrich2...'
csv.field_size_limit(sys.maxsize)
def import_text(filename, separator):
    for line in csv.reader(open(filename), delimiter=separator, 
                           skipinitialspace=True):
        if line:
            yield line
            
################################################################################
# build dictionaries to convert repclass and rep families'
if is_bed == "FALSE":
	repeatclass = {}
	repeatfamily = {}
	fin = import_text(annotation_file, ' ')
	x = 0
	for line in fin:
	    if x>2:
	        classfamily =[]
	        classfamily = line[10].split(os.path.sep)
	        line9 = line[9].replace("(","_").replace(")","_").replace("/","_")
	        repeatclass[line9] = classfamily[0]
	        if len(classfamily) == 2:
	            repeatfamily[line9] = classfamily[1]
	        else:
	            repeatfamily[line9] = classfamily[0]
	    x +=1
if is_bed == "TRUE":
	repeatclass = {}
	repeatfamily = {}
	fin = open(annotation_file, 'r')
	for line in fin:
		line=line.strip('\n')
		line=line.split('\t')
	        theclass =line[4]
	        thefamily = line[5]
	        line3 = line[3].replace("(","_").replace(")","_").replace("/","_")
	        repeatclass[line3] = theclass 
                repeatfamily[line3] = thefamily
fin.close()

################################################################################
# build list of repeats initializing dictionaries for downstream analysis'
fin = import_text(setup_folder + os.path.sep + 'repgenomes_key.txt', '\t')
repeat_key ={}
rev_repeat_key ={}
repeat_list = []
reptotalcounts = {}
classfractionalcounts = {}
familyfractionalcounts = {}
classtotalcounts = {}
familytotalcounts = {}
reptotalcounts_simple = {}
fractionalcounts = {}
i = 0
for line in fin:
    reptotalcounts[line[0]] = 0
    fractionalcounts[line[0]] = 0
    if repeatclass.has_key(line[0]):
	classtotalcounts[repeatclass[line[0]]] = 0
	classfractionalcounts[repeatclass[line[0]]] = 0
    if repeatfamily.has_key(line[0]):
	familytotalcounts[repeatfamily[line[0]]] = 0
	familyfractionalcounts[repeatfamily[line[0]]] = 0
    if repeatfamily.has_key(line[0]):
        if repeatfamily[line[0]] == simple_repeat:
            reptotalcounts_simple[simple_repeat] = 0
    else:
        reptotalcounts_simple[line[0]] = 0
    repeat_list.append(line[0])
    repeat_key[line[0]] = int(line[1])
    rev_repeat_key[int(line[1])] = line[0]
fin.close()
################################################################################
# map the repeats to the psuedogenomes:
if not os.path.exists(outputfolder):
	os.mkdir(outputfolder)
################################################################################
# Conduct the regions sorting
print 'Conducting region sorting on unique mapping reads....'
fileout= outputfolder + os.path.sep + outputfile_prefix + '_regionsorter.txt'
with open(fileout, 'w') as stdout:
	#check version of bedtools to determine which syntax to use for coverageBed
	if (int(v1) == 2 and int(v2) < 24):
		command = shlex.split("coverageBed -abam " +unique_mapper_bam+" -b " +setup_folder + os.path.sep + 'repnames.bed')
	else:
		command = shlex.split("coverageBed -a " + setup_folder + os.path.sep + 'repnames.bed' + " -b " + unique_mapper_bam)
	p = subprocess.Popen(command, stdout=stdout)
p.communicate()
stdout.close()
filein = open(outputfolder + os.path.sep + outputfile_prefix + '_regionsorter.txt','r')
counts = {}
sumofrepeatreads=0
for line in filein:
	line= line.split('\t')
	if not counts.has_key(str(repeat_key[line[3]])):
		counts[str(repeat_key[line[3]])]=0
	counts[str(repeat_key[line[3]])]+=int(line[4])
	sumofrepeatreads+=int(line[4])
print 'Identified ' + str(sumofrepeatreads) + ' unique reads that mapped to repeats.'
################################################################################
if paired_end == 'TRUE':
	if not os.path.exists(outputfolder + os.path.sep + 'pair1_'):
		os.mkdir(outputfolder + os.path.sep + 'pair1_')
	if not os.path.exists(outputfolder + os.path.sep + 'pair2_'):
		os.mkdir(outputfolder + os.path.sep + 'pair2_')
	folder_pair1 = outputfolder + os.path.sep + 'pair1_'
	folder_pair2 = outputfolder + os.path.sep + 'pair2_'
################################################################################
	print "Processing repeat psuedogenomes..."
	ps = []
	psb= []
	ticker= 0
	for metagenome in repeat_list:
   		metagenomepath = setup_folder + os.path.sep + metagenome
		file1=folder_pair1 + os.path.sep + metagenome + '.txt'
		file2 =folder_pair2 + os.path.sep + metagenome + '.txt'
		with open(file1, 'w') as stdout:
			p = subprocess.Popen('bowtie2 ' + b_opt + ' ' + '-x ' + metagenomepath + ' ' + fastqfile_1 + ' | grep \"repname\" -', shell=True, stdout=stdout)
			
		with open(file2, 'w') as stdout:
			pp = subprocess.Popen('bowtie2 ' + b_opt + ' ' + '-x ' + metagenomepath + ' ' + fastqfile_2 + ' | grep \"repname\" -', shell=True, stdout=stdout)
		ps.append(p)
		ticker +=1
		psb.append(pp)
		ticker +=1
		if ticker == cpus:
			for p in ps:
				p.communicate()
			for p in psb:
				p.communicate()
			ticker = 0
			psb =[]
			ps = []
	if len(ps) > 0:
		for p in ps:
			p.communicate()
	stdout.close()
    
################################################################################
# combine the output from both read pairs:
	print 'sorting and combining the output for both read pairs...'
	if not os.path.exists(outputfolder + os.path.sep + 'sorted_'):
		os.mkdir(outputfolder + os.path.sep + 'sorted_')
	sorted_ = outputfolder + os.path.sep + 'sorted_'
	for metagenome in repeat_list:
		file1 = folder_pair1 + os.path.sep + metagenome + '.txt'
		file2 = folder_pair2 + os.path.sep + metagenome + '.txt'
		fileout= sorted_ + os.path.sep + metagenome + '.txt'
		with open(fileout, 'w') as stdout:
			p1 = subprocess.Popen(['cat',file1,file2], stdout = subprocess.PIPE)
			p2 = subprocess.Popen(['cut', '-f1',"-d "], stdin = p1.stdout, stdout = subprocess.PIPE)
			p3 = subprocess.Popen(['cut', '-f1', "-d/"], stdin = p2.stdout, stdout = subprocess.PIPE)
			p4 = subprocess.Popen(['sort'], stdin=p3.stdout, stdout = subprocess.PIPE)
			p5 = subprocess.Popen(['uniq'], stdin=p4.stdout, stdout = stdout)
			p5.communicate()
		stdout.close()
	print 'completed ...'
################################################################################
if paired_end == 'FALSE':
	if not os.path.exists(outputfolder + os.path.sep + 'pair1_'):
		os.mkdir(outputfolder + os.path.sep + 'pair1_')
	folder_pair1 = outputfolder + os.path.sep + 'pair1_'
################################################################################
	ps = []
	ticker= 0
	print "Processing repeat psuedogenomes..."
	for metagenome in repeat_list:
    		metagenomepath = setup_folder + os.path.sep + metagenome
		file1=folder_pair1 + os.path.sep + metagenome + '.txt'
		with open(file1, 'w') as stdout:
			p = subprocess.Popen('bowtie2 ' + b_opt + ' ' + '-x ' + metagenomepath + ' ' + fastqfile_1 + ' | grep \"repname\" -', shell=True, stdout=stdout)
		ps.append(p)
		ticker +=1
		if ticker == cpus:
			for p in ps:
				p.communicate()
			ticker = 0
			ps = []
	if len(ps) > 0:
		for p in ps:
			p.communicate()
	stdout.close()
    
################################################################################
# combine the output from both read pairs:
	print 'Sorting and combining the output for both read pairs....'
	if not os.path.exists(outputfolder + os.path.sep + 'sorted_'):
		os.mkdir(outputfolder + os.path.sep + 'sorted_')
	sorted_ = outputfolder + os.path.sep + 'sorted_'
	for metagenome in repeat_list:
		file1 = folder_pair1 + os.path.sep + metagenome + '.txt'
		fileout= sorted_ + os.path.sep + metagenome + '.txt'
		with open(fileout, 'w') as stdout:
			p1 = subprocess.Popen(['cat',file1], stdout = subprocess.PIPE)
			p2 = subprocess.Popen(['cut', '-f1'], stdin = p1.stdout, stdout = subprocess.PIPE)
			p3 = subprocess.Popen(['cut', '-f1', "-d/"], stdin = p2.stdout, stdout = subprocess.PIPE)
			p4 = subprocess.Popen(['sort'], stdin = p3.stdout,stdout = subprocess.PIPE)
			p5 = subprocess.Popen(['uniq'], stdin = p4.stdout,stdout = stdout)
			p5.communicate()
		stdout.close()
	print 'completed ...'
    
################################################################################
# build a file of repeat keys for all reads
print 'Writing and processing intermediate files...'
sorted_ = outputfolder + os.path.sep + 'sorted_'
readid = {}
sumofrepeatreads=0
for rep in repeat_list: 
    for data in import_text(sorted_ + os.path.sep + rep + '.txt', '\t'):
	readid[data[0]] = ''
for rep in repeat_list: 
    for data in import_text(sorted_ + os.path.sep + rep + '.txt', '\t'):	
	readid[data[0]]+=str(repeat_key[rep]) + str(',')
for subfamilies in readid.values():
    if not counts.has_key(subfamilies):
	counts[subfamilies]=0
    counts[subfamilies] +=1
    sumofrepeatreads+=1
del readid
print 'Identified ' + str(sumofrepeatreads) + ' reads that mapped to repeats for unique and multimappers.'

################################################################################
print "Conducting final calculations..."
# build a converter to numeric label for repeat and yield a combined list of repnames seperated by backslash
def convert(x):
    x = x.strip(',')
    x = x.split(',')
    global repname
    repname = ""
    for i in x:
        repname = repname + os.path.sep + rev_repeat_key[int(i)] 
# building the total counts for repeat element enrichment...
for x in counts.keys():
    count= counts[x]
    x = x.strip(',')
    x = x.split(',')
    for i in x:
        reptotalcounts[rev_repeat_key[int(i)]] += int(count)
# building the fractional counts for repeat element enrichment...
for x in counts.keys():
    count= counts[x]
    x = x.strip(',')
    x = x.split(',')
    splits = len(x)
    for i in x:
        fractionalcounts[rev_repeat_key[int(i)]] += float(numpy.divide(float(count),float(splits)))
# building categorized table of repeat element enrichment... 
repcounts = {}
repcounts['other'] = 0
for key in counts.keys():
        convert(key)
        repcounts[repname] = counts[key]
# building the total counts for class enrichment...
for key in reptotalcounts.keys():
	classtotalcounts[repeatclass[key]] += reptotalcounts[key]
# building total counts for family enrichment...
for key in reptotalcounts.keys():
	familytotalcounts[repeatfamily[key]] += reptotalcounts[key]
# building unique counts table'
repcounts2 = {}
for rep in repeat_list:
    if repcounts.has_key("/" +rep):
        repcounts2[rep] = repcounts["/" +rep]
    else:
        repcounts2[rep] = 0
# building the fractionalcounts counts for class enrichment...
for key in fractionalcounts.keys():
	classfractionalcounts[repeatclass[key]] += fractionalcounts[key]
# building fractional counts for family enrichment...
for key in fractionalcounts.keys():
	familyfractionalcounts[repeatfamily[key]] += fractionalcounts[key]   
    
################################################################################
print 'Writing final output and removing intermediate files...'
# print output to file of the categorized counts and total overlapping counts:
if allcountmethod  == "TRUE":
	fout1 = open(outputfolder + os.path.sep + outputfile_prefix + '_total_counts.txt' , 'w')
	for key in reptotalcounts.keys():
	    print >> fout1, str(key) + '\t' + repeatclass[key] + '\t' + repeatfamily[key] + '\t' + str(reptotalcounts[key])
	fout2 = open(outputfolder + os.path.sep + outputfile_prefix + '_class_total_counts.txt' , 'w')
	for key in classtotalcounts.keys():
	    print >> fout2, str(key) + '\t' + str(classtotalcounts[key])  
	fout3 = open(outputfolder + os.path.sep + outputfile_prefix + '_family_total_counts.txt' , 'w')
	for key in familytotalcounts.keys():
	    print >> fout3, str(key) + '\t' + str(familytotalcounts[key])  
	fout4 = open(outputfolder + os.path.sep + outputfile_prefix + '_unique_counts.txt' , 'w')
	for key in repcounts2.keys():
	    print >> fout4, str(key) + '\t' + repeatclass[key] + '\t' + repeatfamily[key] + '\t' + str(repcounts2[key])     
    	fout5 = open(outputfolder + os.path.sep + outputfile_prefix + '_class_fraction_counts.txt' , 'w')
	for key in classfractionalcounts.keys():
	    print >> fout5, str(key) + '\t' + str(classfractionalcounts[key])  
	fout6 = open(outputfolder + os.path.sep + outputfile_prefix + '_family_fraction_counts.txt' , 'w')
	for key in familyfractionalcounts.keys():
	    print >> fout6, str(key) + '\t' + str(familyfractionalcounts[key])
	fout7 = open(outputfolder + os.path.sep + outputfile_prefix + '_fraction_counts.txt' , 'w')
	for key in fractionalcounts.keys():
	    print >> fout7, str(key) + '\t' + repeatclass[key] + '\t' + repeatfamily[key] + '\t' + str(int(fractionalcounts[key]))
    	fout1.close()
	fout2.close()
	fout3.close()
	fout4.close()
	fout5.close()
	fout6.close()
	fout7.close()
else:
	fout1 = open(outputfolder + os.path.sep + outputfile_prefix + '_class_fraction_counts.txt' , 'w')
	for key in classfractionalcounts.keys():
	    print >> fout1, str(key) + '\t' + str(classfractionalcounts[key])  
	fout2 = open(outputfolder + os.path.sep + outputfile_prefix + '_family_fraction_counts.txt' , 'w')
	for key in familyfractionalcounts.keys():
	    print >> fout2, str(key) + '\t' + str(familyfractionalcounts[key])
	fout3 = open(outputfolder + os.path.sep + outputfile_prefix + '_fraction_counts.txt' , 'w')
	for key in fractionalcounts.keys():
	    print >> fout3, str(key) + '\t' + repeatclass[key] + '\t' + repeatfamily[key] + '\t' + str(int(fractionalcounts[key]))
	fout1.close()
	fout2.close()
	fout3.close()
    
################################################################################
# Remove Large intermediate files
if os.path.exists(outputfolder + os.path.sep + outputfile_prefix + '_regionsorter.txt'):
	os.remove(outputfolder + os.path.sep + outputfile_prefix + '_regionsorter.txt')
if os.path.exists(outputfolder + os.path.sep + 'pair1_'):
	shutil.rmtree(outputfolder + os.path.sep + 'pair1_')
if os.path.exists(outputfolder + os.path.sep + 'pair2_'):
	shutil.rmtree(outputfolder + os.path.sep + 'pair2_')
if os.path.exists(outputfolder + os.path.sep + 'sorted_'):
	shutil.rmtree(outputfolder + os.path.sep + 'sorted_')

print "... Done"