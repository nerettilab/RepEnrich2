#RepEnrich2
## Tutorial By Nicholas Skvir  
Email: [nicholas_skvir@brown.edu](mailto:nicholas_skvir@brown.edu)  

### Dependencies
This example is for mouse genome **mm9**. Before getting started you should make sure you have installed the dependencies
for RepEnrich2.  

RepEnrich2 currently requires:  
[Python](https://www.python.org/downloads/) version 2.7+,  
[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),  
[Bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html),  
[Samtools](http://www.htslib.org/),  
[BioPython](http://biopython.org)  

BioPython is used by the python scripts, and can be installed with the following command:

    pip install BioPython

Currently I am using Python version 2.7.12+, Bedtools version 2.25.0, Bowtie2 version 2.2.9, Samtools version 1.3.1, and BioPython version 1.66-py2.7.  

Lastly, RepEnrich2 requires a Bowtie2 indexed genome in fasta format available. (Example `mm9.fa`)  


### Step 1) Attain repetitive element annotation
I have temporarily provided the setup for the human genome (build hg19 and hg38) and the mouse genome (build mm9) available [here] (https://drive.google.com/folderview?id=0B1dD8MQRH4qZfmlxOGwtSXRnWDFaVldqbkExdXItZGpySm1mVmhlTVladThHWWhGMmxrLTQ&usp=sharing). After downloading you can extract the files using:

    gunzip hg19_repeatmasker_clean.txt.gz
    tar -zxvf RepEnrich_setup_hg19.tar.gz
  
To yield hg19_repeatmasker_clean.txt annotation file and RepEnrich_setup_hg19 setup folder.  The annotation files I am using are repeatmasker files with simple and low-complexity repeats removed (satellite repeats and transposons are still present).  If you choose to use these files for the set-up you can skip ahead to step 2.

The RepEnrich2_setup script will build the annotation required by RepEnrich. The default is a repeatmasker file which can be
downloaded from [repeatmasker.org](http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html),
(for instance, find the `mm9.fa.out.gz` download
[here](http://www.repeatmasker.org/genomes/mm9/RepeatMasker-rm328-db20090604/mm9.fa.out.gz).
Once you have downloaded the file you can unzip it and rename it:

    gunzip mm9.fa.out.gz
    mv mm9.fa.out mm9_repeatmasker.txt

This is what the file looks like:

	SW perc perc perc query position in query matching repeat position in repeat score div. del. ins. sequence begin end (left) repeat class/family begin end  (left) ID
	687 17.4 0.0 0.0 chr1 3000002 3000156 (194195276) C L1_Mur2 LINE/L1 (4310) 1567 1413 1
	917 21.4 11.4 4.5 chr1 3000238 3000733 (194194699) C L1_Mur2 LINE/L1 (4488) 1389 913 1
	845 23.3 7.6 11.4 chr1 3000767 3000792 (194194640) C L1_Mur2 LINE/L1 (6816) 912 887 1
	621 25.0 6.5 3.7 chr1 3001288 3001583 (194193849) C Lx9 LINE/L1 (1596) 6048 5742 3

The RepEnrich2_setup script will also allow you to build the annotation required by RepEnrich2 for a custom set of elements using a bed file. So if you want to examine mm9 LTR repetitive elements; you can build this file using the the repeatmasker track from [UCSC genome table browser](http://genome.ucsc.edu/cgi-bin/hgTables).

To do this, select genome `mm9`, click the edit box next to _Filter_, fill in the repclass does match with `LTR`, then click submit. Back at the table browser select option `Selected fields from primary and related tables`, name the output file something like `mm9_LTR_repeatmasker.bed`, and click `Get output`. On the next page select `genoName`, `genoStart`, `genoEnd`, `repName`, `repClass`, `repFamily` then download the file.

The UCSC puts a header on the file that needs to be removed: 

    tail -n +3 mm9_LTR_repeatmasker.bed | head -n -4 > mm9_LTR_repeatmasker_fix.bed
    mv mm9_LTR_repeatmasker_fix.bed mm9_LTR_repeatmasker.bed

This is what our custom mm9 LTR retrotransposon bed file looks like:

	$ head mm9_LTR_repeatmasker.bed

	chr1 3001722 3002005 RLTR25A LTR ERVK
	chr1 3002051 3002615 RLTR25A LTR ERVK
	chr1 3016886 3017193 RLTRETN_Mm LTR ERVK
	chr1 3018338 3018653 RLTR14 LTR ERV1

Note: It is important to get the column format right: 

* Column 1: Chromosome
* Column 2: Start
* Column 3: End
* Column 4: Repeat_name
* Column 5: Class
* Column 6: Family

The file should be tab delimited. If there is no information on class or family, you can replace these columns with the repeat name or an arbitrary label such as `group1`.


### Step 2) Run the setup for RepEnrich2

Now that we have our annotation files we can move on to running the setup for RepEnrich2. First load the dependencies (if you use Environment Modules - otherwise just make sure that these programs are available in your `PATH`).

    module load bowtie2
    module load bedtools
    module load samtools

Next run the setup using the type of annotation you have selected (default):

    python RepEnrich2_setup.py /data/mm9_repeatmasker.txt /data/mm9.fa /data/setup_folder_mm9

custom bed file:

    python RepEnrich2_setup.py /data/mm9_LTR_repeatmasker.bed /data/mm9.fa /data/setup_folder_mm9 --is_bed TRUE

The previous commands have setup RepEnrich2 annotation that is used in downstream analysis of data. You only have to do the setup step once for an organism of interest. One cautionary note is that RepEnrich2 is only as reliable as the genome annotation of repetitive elements for your organism of interest. Therefore, RepEnrich2 performance may not be optimal for poorly annotated genomes. 


### Step 3) Map the data to the genome

After the setup of the RepEnrich2 we now have to map our data uniquely to the genome before running. This is because RepEnrich2 treats uniquely mapping and multi-mapping reads separately. This requires first mapping the data to the genome using Bowtie2, then running the RepEnrich2_subset.py script to segregate the reads into separate files. 

First, perform mapping using Bowtie2 (example using paired end data):

     bowtie2 -q -p 16 -x /path/to/annotation/mm9/mm9 -1 /home/Dataset/Sample_L001_R1.fastq -2 /home/Dataset/Sample_L001_R2.fastq -S /home/output/sample_mapped.sam

An explanation of the Bowtie2 options (detailed explanations found [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#command-line)):
* `bowtie2 `
* `-q` -- reads are fastq files (extension .fq or .fastq)
* `-p 16` -- using 16 cpus
* `-x <bt2-index>` -- Basename of index for reference genome. (Name of the index files up to but not including the .1.bt2 / etc.)
* `-1 <file_path_1>` -- path to first paired end fastq file (use -U instead of -1/-2 if using unpaired data)
* `-2 <file_path_2>` -- path to second paired end fastq file
* `-S <output_file> ` --output a Sam file with the specified filename/path

The `.out` file from this run should inform you of the total mapping reads for your alignment



The .sam file should then be converted to a .bam file using Samtools:

    samtools view -bS sample.sam > sample.bam

Lastly, the RepEnrich2_subset.py script should be run to output discrete files for uniquely and multi-mapping reads:

    python RepEnrich2_subset.py /path/to/sample.bam 30 Sample_Name --pairedend TRUE

An explanation of the RepEnrich2_subset.py commands:
* `<bamfile_to_subset>` -- path to the bamfile generated in the previous step
* `30` -- the mapq value on which to subset uniquely mapping reads from multi-mapping reads; 30 is recommended for bowtie2
* `<sample_name>` -- a prefix for the output files
* `--pairedend TRUE` -- specifies whether the data is paired end, default is FALSE if unspecified



The script should output a .bam file of uniquely mapping reads and .fastq file(s) of multimapping reads (sample_name_unique.bam, sample_name_multimap_R1.fastq, sample_name_multimap_R2.fastq for example).



### Step 4) Run RepEnrich2 on the data

Now we have all the information we need to run RepEnrich2.
Here is an example (for default annotation):

    python RepEnrich2.py /data/mm9_repeatmasker.txt /data/Sample_Output_Folder/ Sample_Name /data/RepEnrich2_setup_mm9 /data/sample_name_multimap_R1.fastq --fastqfile2 /data/sample_name_multimap_R2.fastq /data/sample_name_unique.bam --cpus 16 --pairedend TRUE

for custom bed file annotation:

    python RepEnrich2.py /data/mm9_LTR_repeatmasker.bed /data/Sample_Output_Folder Sample_Name /data/RepEnrich2_setup_mm9 /data/sample_name_multimap_R1.fastq --fastqfile2 /data/sample_name_multimap_R2.fastq /data/sample_name_unique.bam --is_bed TRUE --cpus 16 --pairedend TRUE

An explanation of the RepEnrich2 commands:

	python RepEnrich2.py 
		<repeat_annotation>
		<output_folder> -- location to output files; creates a new folder if one doesn't exist
		<output_prefix>
		<RepEnrich2_setup_folder>
		<multimapping_read_1.fastq>
    (--fastqfile2) <multimapping_read_2.fastq>
		<unique_mapping_reads.bam>
		(--is_bed TRUE) -- if using a custom .bed file for repeat annotation
		(--cpus 16)
    (--pairedend TRUE) -- if using paired end data


### Step 5) Processing the output of RepEnrich

The final outputs will be in the path `/data/Sample_Name`. This will include several files. The most important of which is the `sampleA_fraction_counts.txt` file. This is the estimated counts for the repeats. I use this file to build a table of counts for all my conditions (by pasting the individual `*_fraction_counts.txt` files together for my complete experiment).

You can use the compiled counts file to do differential expression analysis similar to what is done for genes. We use [EdgeR](http://www.bioconductor.org/packages/release/bioc/html/edgeR.html) or [DESEQ](http://bioconductor.org/packages/release/bioc/html/DESeq.html) to do the differential expression analysis. These are R packages that you can download from [bioconductor](http://bioconductor.org/).

When running the EdgeR differential expression analysis you can follow the examples in the EdgeR manual. I manually input the
library sizes (the total mapping reads we obtained in the tutorial). Some of the downstream analysis, though, is left to your discretion. There are multiple ways you can do the differential expression analysis. I use the `GLM` method within the EdgeR packgage, although DESeq has similar methods and EdgeR also has a more straightforward approach called `exactTest`. Below is a sample EdgeR script used to do the differential analysis of repeats for young, old, and very old mice. The file `counts.csv` contains the ouput from RepEnrich that was made by pasting the individual `*_fraction_counts.txt` files together for my complete experiment. 


## Example Script for EdgeR differential enrichment analysis

```r
# EdgeR example

# Setup - Install and load edgeR
source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
library('edgeR')

# In the case of a pre-assembled file of the fraction count output do the following:
# counts <- read.csv(file = "counts.csv")

# In the case of seperate outputs, load the RepEnrich results - fraction counts
young_r1 <- read.delim('young_r1_fraction_counts.txt', header=FALSE)
young_r2 <- read.delim('young_r2_fraction_counts.txt', header=FALSE)
young_r3 <- read.delim('young_r3_fraction_counts.txt', header=FALSE)
old_r1 <- read.delim('old_r1_fraction_counts.txt', header=FALSE)
old_r2 <- read.delim('old_r2_fraction_counts.txt', header=FALSE)
old_r3 <- read.delim('old_r3_fraction_counts.txt', header=FALSE)
v_old_r1 <- read.delim('veryold_r1_fraction_counts.txt', header=FALSE)
v_old_r2 <- read.delim('veryold_r2_fraction_counts.txt', header=FALSE)
v_old_r3 <- read.delim('veryold_r3_fraction_counts.txt', header=FALSE)

#' Build a counts table
counts <- data.frame(
  row.names = young_r1[,1],
  young_r1 = young_r1[,4], young_r2 = young_r2[,4], young_r3 = young_r3[,4],
  old_r1 = old_r1[,4], old_r2 = old_r2[,4], old_r3 = old_r3[,4],
  v_old_r1 = v_old_r1[,4], v_old_r2 = v_old_r2[,4], v_old_r3 = v_old_r3[,4]
)

# Build a meta data object. I am comparing young, old, and veryold mice.
# I manually input the total mapping reads for each sample.
# The total mapping reads are calculated using the bowtie logs:
# # of reads processed - # reads that failed to align
meta <- data.frame(
	row.names=colnames(counts),
	condition=c("young","young","young","old","old","old","veryold","veryold","veryold"),
	libsize=c(24923593,28340805,21743712,16385707,26573335,28131649,34751164,37371774,28236419)
)

# Define the library size and conditions for the GLM
libsize <- meta$libsize
condition <- factor(meta$condition)
design <- model.matrix(~0+condition)
colnames(design) <- levels(meta$condition)

# Build a DGE object for the GLM
y <- DGEList(counts=counts, lib.size=libsize)

# Normalize the data
y <- calcNormFactors(y)
y$samples
plotMDS(y)

# Estimate the variance
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y)

# Build an object to contain the normalized read abundance
logcpm <- cpm(y, log=TRUE, lib.size=libsize)
logcpm <- as.data.frame(logcpm)
colnames(logcpm) <- factor(meta$condition)

# Conduct fitting of the GLM
yfit <- glmFit(y, design)

# Initialize result matrices to contain the results of the GLM
results <- matrix(nrow=dim(counts)[1],ncol=0)
logfc <- matrix(nrow=dim(counts)[1],ncol=0)

# Make the comparisons for the GLM
my.contrasts <- makeContrasts(
	veryold_old = veryold – old,
	veryold_young = veryold – young,
	old_young = old – young,
	levels = design
)

# Define the contrasts used in the comparisons
allcontrasts = c(
	"veryold_old",
	"veryold_young",
	"old_young"
)

# Conduct a for loop that will do the fitting of the GLM for each comparison
# Put the results into the results objects
for(current_contrast in allcontrasts) {
	lrt <- glmLRT(yfit, contrast=my.contrasts[,current_contrast])
	plotSmear(lrt, de.tags=rownames(y))
	title(current_contrast)
	res <- topTags(lrt,n=dim(c)[1],sort.by="none")$table
	colnames(res) <- paste(colnames(res),current_contrast,sep=".")
	results <- cbind(results,res[,c(1,5)])
	logfc <- cbind(logfc,res[c(1)])
}

# Add the repeat types back into the results.
# We should still have the same order as the input data
results$class <- young_r1[,2]
results$type <- young_r1[,3]

# Sort the results table by the logFC
results <- results[with(results, order(-abs(logFC.old_young))), ]

# Save the results
write.table(results, 'results.txt', quote=FALSE, sep="\t")

# Plot Fold Changes for repeat classes and types
for(current_contrast in allcontrasts) {
  logFC <- results[, paste0("logFC.", current_contrast)]
  # Plot the repeat classes
  classes <- with(results, reorder(class, -logFC, median))
  par(mar=c(6,10,4,1))
  boxplot(logFC ~ classes, data=results, outline=FALSE, horizontal=TRUE,
          las=2, xlab="log(Fold Change)", main=current_contrast)
  abline(v=0)
  # Plot the repeat types
  types <- with(results, reorder(type, -logFC, median))
  boxplot(logFC ~ types, data=results, outline=FALSE, horizontal=TRUE,
          las=2, xlab="log(Fold Change)", main=current_contrast)
  abline(v=0)
}

```


Note that the objects `logfc` contains the differential expression for the contrast, `logcpm` contains the normalized read abundance, and `result` contains both the differential expression and the false discovery rate for the experimental comparison. I recommended reading more about these in the [EdgeR manual](http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf).

