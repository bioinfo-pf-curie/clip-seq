# Clip Quick Start Guide

See NEWS for information about changes in this and previous versions.

## What is the CLIP-seq pipeline ?

CLIP-seq was set up to process CLIP sequencing data from aligned sequencing reads. It is a pipeline including two parts: <br /> * the reads processing (cleaning, mapping, peak calling and annotation of peaks) <br /> * the differential analysis. <br /> The normalization step in the differential analysis can be defined as the normalization of the number of reads in a peak (in CLIP-seq samples) by the number of reads in the gene (in RNA-seq sample) where the peak is located on the gene.

## Contact

For any questions about the pipeline, please contact <mandy.cadix@curie.fr>

## Requirement

The following dependancies are required :

* R (>3.2.0) with the *RColorBrewer*, *ggplot2*, *rtracklayer*, *DESeq2*, *magrittr*, *dplyr*, *qplots*, *Rsubread*, *plyr* and *GenomicRanges* packages
* Samtools (>1.1)
* BEDTools (>2.21.0)
* Python (>2.7.9) with the *os* and *argparse* packages
* Cutadapt (>1.8.2)
* Tophat (>2.0.14)
* Java/jkd (>1.7.0_45)
* Picard_tools
* Piranha (>1.2.1)

## How to install it ?

To install the CLIP-seq pipeline, simply extract the archive and set up the configuration file with the paths to dependencies.

    TO DO

## Config File

* Annotation Files

The pipeline is using a annotation files with gene annotation information. These file is based on UCSC Refseq gene. WINDOW\_SIZE is a parameter to select a part of around the end of gene (upstream or downstream at the end of gene), this size is pb unit. In order to generate all required annotation files, please set the ANNOTATION\_DIR, ORG, UCSC\_EXPORT and WINDOW\_SIZE in the configuration file.

    BUILD_ANNOT=1
    ORG=hg19
    UCSC_EXPORT=refseq_export_hg19.csv
    WINDOW_SIZE=100


* Trimming

To treat reads with adapter, it is necessary to know side adapters 5' and 3':

    ADAPTER5="GTTCAGAGTTCTACAGTCCGACGATC"
    ADAPTER3="TGGAATTCTCGGGTGCCAAGG"


## How to use it ?

* PART 1 : reads processing

CLIP-seq can be used for a single sample. In order to use the pipeline, please set up the configuration according to your analysis, and run the following command to do the reads processing:

    ./script/pipeline_clip-seq.bash -c CONFIG -i INPUT_FILE -s STEP -n SAMPLE_NUMBER -o OUTPUT_DIR


**VARIABLE NAME** | **CONTENT**
----------------- | -----------
-c   | The TXT format configuration file
-i   | The fastq format file of the CLIP-seq sample
-s   | The steps of the pipeline that are to be launched
-n   | The sample identifiers (only numbers)
-o   | The output directory
-h   | Print help message


There are 15 differents steps (-s): 

**VARIABLE NAME**             | **CONTENT**
-----------------             | -----------
clean\_fastq                  | Remove reads which have had problems with the sequencer. CAUTION: only the files in fastq format are accepted
trimming                      | Remove adapters and to select reads with a minimum of nucleotides
mapping                       | Map reads on tophat2
MAPQ                          | Use the mapping quality to select reads. On this step we create a sort file
duplicate\_reads              | Remove (or only mark) the duplicated reads. On this step we create a index file
peak\_calling                 | Detect of peaks
annotation\_gene              | Build an annotation file
treatment\_peaks\_file        | Adapt the peak file for the annotation step
annotation\_peaks             | Annotate peaks
annotation\_peaks\_last\_exon | Annotate peaks present in the last exon region
annotation\_peaks\_downstream | Annotate peaks present in the downstream region from the end of the gene
annotation\_peaks\_upstream   | Annotate peaks present in the upstream region from the end of the gene
annotation\_peaks\_extended   | Annotate peaks present in all gene and the downstream region from the end of the gene
all                           | Launch all previous step


* PART 2 : differential analysis

If you want to compare treated VS untreated please enter treated samples before untreated samples (for RNA-seq, CLIP-seq and group) and run the following command to do the differential analysis:


    ./script/pipeline_differential_analysis.bash -c CONFIG -l INPUT_LIST -b BAM_RNA_FILES -s STEP -n SAMPLE_TO_COMBINE -g SAMPLE_GROUP -o OUTPUT_DIR

**VARIABLE NAME** | **CONTENT**
----------------- | -----------
 -c   | The TXT format configuration file
 -l   | The TXT format file of BED files, obtained from the previous script: <br /> id\_CLIP\_sample\<TAB\>path\_of\_bed\_file\<TAB\>condition\_of\_this\_sample <br /> condition\_of\_this\_sample is the condition of the sample (e.g: Untreated, treated ...) CAUTION: Condition must be write with letters (no number or symbol)
 -b   | The TXT format file of BAM file sorted of RNA-seq samples: <br /> id\_CLIP\_sample\<TAB\>path\_of\_bam\_file\<TAB\>condition\_of\_this\_sample
 -s   | The steps of the pipeline that are to be launched
 -n   | The identifier of CLIP-seq samples (numbers) separated by comma
 -g   | The identifier of sample condition separated by comma (0 for untreated sample and 1 for treated sample)
 -o   | The output directory
 -h   | Print help message


There are 7 differents steps (-s): 

**VARIABLE NAME**          | **CONTENT**
-----------------          | -----------
concatenate\_peaks         | Concatenate peaks present in different CLIP-seq files
merge\_peaks               | Merge peaks between differents CLIP-seq samples
tableOfCounts\_peaks       | Create a table of counts of peaks detected in CLIP-seq samples
tableOfCounts\_genes       | Create a table of counts of genes detected in RNA-seq samples
concatenate\_peaks\_genes  | Concatenate tables of counts of CLIP-seq and RNA-seq in the same file
differential\_analysis     | Make the differential analysis
color\_peak                | Create a file with rgb code color per peak
all                        | Launch all previous step


## Test dataset

    TO DO

