CLIP-seq Quick Start Guide
==========================

See NEWS for information about changes in this and previous versions

What is the CLIP-seq pipeline ?
-------------------------------

CLIP-seq was set up to process CLIP sequencing data from sequencing reads. This pipeline is developed for **single-end** data.
It includes the trimming, mapping, removing duplicated reads, peak identification, annotation dans differential analysis.

Contact
-------

For any questions about the pipeline, please contact <mandy.cadix@curie.fr>

How to install it?
------------------

The following dependancies are required :

* [R] (https://www.r-project.org/) (version 3.4.0) with the *RColorBrewer (v1.1-2)*, *ggplot2 (v2.2.1)*, *rtracklayer (v1.36.3)*, *DESeq2 (v1.16.1)*, *magrittr (v1.5)*, *dplyr (v0.7.2)*, *gplots (v3.0.1)*, *plyr (v1.8.4)*, *Rsubread (v1.28.0)* and *GenomicRanges (v1.28.3)* packages

* [Samtools] (http://samtools.sourceforge.net) (version 1.1)

* [BEDTools] (http://bedtools.readthedocs.io/en/latest/content/installation.html) (version 2.25.0)

* [Python] (https://www.python.org/downloads/release/python-279) (version 2.7.12) with the *os*, *argparse (v1.1)*, *re (v2.2.1)* and *gzip* packages

* [Cutadapt] (http://cutadapt.readthedocs.io/en/stable/installation.html) (version 1.12)

* [Tophat2] (http://ccb.jhu.edu/software/tophat/downloads/) (version 2.1.1)

* [Bowtie2] (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) (version 2.2.6)

* [Piranha] (https://github.com/smithlabcode/piranha) (version 1.2.1)

* [Java/jdk] (http://www.oracle.com/technetwork/java/javase/downloads/index-jsp-138363.html) (version 1.7.0\_45)

* [Fastqc] (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.5)

* [Picard\_tools] (https://broadinstitute.github.io/picard/) (version 1.97)


Annotation Files
----------------

The pipeline is using a annotation files with gene annotation information. These file is based on UCSC Refseq gene. WINDOW\_SIZE is a parameter to select a part of around the end of gene (upstream or downstream at the end of gene), this size is bp unit. In order to generate all required annotation files, please set the ANNOTATION\_DIR, ORG, UCSC\_EXPORT and WINDOW\_SIZE in the configuration file.

    BUILD_ANNOT=1
    ORG=hg19
    UCSC_EXPORT=refseq_export_hg19.csv
    WINDOW_SIZE=100


Trimming
--------

To treat reads with adapter, it is necessary to know side adapters 5' and 3':

    ADAPTER5="GTTCAGAGTTCTACAGTCCGACGATC"
    ADAPTER3="TGGAATTCTCGGGTGCCAAGG"


How to use it ?
---------------

**PART 1 : reads processing** CLIP-seq can be used for a single sample. In order to use the pipeline, please set up the configuration according to your analysis, and run the following command to do the reads processing:

    ./script/pipeline_clip-seq.bash -c CONFIG -i INPUT_FILE -s STEP -n SAMPLE -o OUTPUT_DIR

    -c: The configuration file
    
    -i: The fastq format file of the CLIP-seq sample. CAUTION: only the files in fastq.gz format are accepted
    
    -s: Steps of this pipeline:
        trimming: Remove reads which have had problems with the sequencer and remove adapters
        fastqc: Quality control for fastq file
        mapping: Map and select reads with tophat2 and the mapping quality.
        rm_dup: Remove or only mark the duplicated reads. On this step we create a sort and index file
        peak_calling: Detect peaks
        annotation: Annotate genes and peaks
        annot_spe: Annotate peaks with a specific window (last exon, downstream & upstream region of peaks, peaks with downstream region)
        all: Launch all previous steps
        
        -n: The sample identifiers (only letter or number; NO symbols)
        
        -o: The output directory
        
        -h: Help
        
        -v: Version


**PART 2 : differential analysis** The normalization step in the differential analysis can be defined as the normalization of the number of reads in a peak (in CLIP-seq samples) by the number of reads in the gene (in RNA-seq sample) where the peak is located on the gene.

If you want to compare treated VS untreated please enter treated samples before untreated samples (for RNA-seq, CLIP-seq and group) and run the following command to do the differential analysis:


    ./script/pipeline_differential_analysis.bash -c CONFIG -l INPUT_LIST -b BAM_RNA_FILES -s STEP -o OUTPUT_DIR

    -c: The configuration file
    
    -l: The list of BED files obtained from the PART 1 (See input_bed_list.txt)
    
    -b: The list of BAM sorted of RNA-seq samples (See input_bam_list.txt)
    
    -s: Steps of this pipeline:
        merge_peaks: Concatenate, sort and merge all peaks present in different samples
        toc_peaks: Create a table of counts for CLIP experiment
        toc_genes: Create a table of counts for RNA experiment
        cat_peaks_genes: Concatenation of tables of counts (CLIP and RNA)
        differential_analysis: Make a differential analysis
        color_peaks: Create color information file. Up and Down regulated peaks (padj <0.05) are colored on red and blue respectively
        all: Launch all previous steps
    
    -o: The output directory

    -h: Help
    
    -v: Version

This is a input list of BED files obtained thanks to the PART 1, for the PART 2, *input_bed_list.txt*

    Sample2  /PATH/SAMPLE2/SAMPLE2_intersectBED_500_bp_DownstreamEnd.bed  treated 1
    Sample4  /PATH/SAMPLE4/SAMPLE4_intersectBED_500_bp_DownstreamEnd.bed  treated 1
    Sample1  /PATH/SAMPLE1/SAMPLE1_intersectBED_500_bp_DownstreamEnd.bed  untreated 0
    Sample3  /PATH/SAMPLE3/SAMPLE3_intersectBED_500_bp_DownstreamEnd.bed  untreated 0
    
This is a input list of BAM RNA-seq files, for the PART 2, *input_bam_list.txt*

    Sample2  /PATH/BAM_RNA-seq_FILE/SAMPLE2_MAPQ_sort.bam treated
    Sample4  /PATH/BAM_RNA-seq_FILE/SAMPLE4_MAPQ_sort.bam treated
    Sample1  /PATH/BAM_RNA-seq_FILE/SAMPLE1_MAPQ_sort.bam untreated
    Sample3  /PATH/BAM_RNA-seq_FILE/SAMPLE3_MAPQ_sort.bam untreated