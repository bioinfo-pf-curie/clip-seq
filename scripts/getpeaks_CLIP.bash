#!/bin/bash

set -o pipefail  # trace ERR through pipes
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value

SOFT="CLIP-seq"
VERSION="2.0.1"

################
##  Function  ##

function evalecho {
    echo 
#    echo $1
    echo
    eval $1
}
################


##  Initialisation  ##
dir=$(dirname $0)

NORMAL="\\033[0;39m"
RED="\\033[1;31m"

die() {
    echo -e "$RED""$*""$NORMAL" 1>&2
    exit 1
}

##  Get arguments  ##
usage()
{
    echo "usage: $0 -c CONF_FILE -i INPUT_FILE -s STEP -n SAMPLE_NUMBER -o OUTPUT_DIR"
    echo "Use option -h for more information"
    exit
}

help()
{
    echo " ***$SOFT  $VERSION *** "
    echo 
    echo "OPTIONS"
    echo "    -c CONFIG : configuration file for $SOFT pipeline"
    echo "    -i INPUT : input file; fastq.gz format"
    echo "    -s STEP : run all or only a subset of the $SOFT workflow"
    echo "	  trimming : Remove reads which have a problem with sequencer and remove adapters"
    echo " 	  fastqc : Quality control for fastq file"
    echo " 	  mapping : Map reads on the reference genome with Tophat2"
    echo " 	  rm_dup : Remove or only mark the duplicated reads"
    echo " 	  peak_calling : Detect peaks"
    echo " 	  annotation : Annotate genes and peaks"
    echo " 	  annot_spe : Annotate peaks with a specific window (last exon, downstream & upstream region of peaks, peaks with downstream region)"
    echo " 	  all : Launch all previous steps"
    echo "    -n SAMPLE : Sample id (letter or letter and number; NO symbols)"
    echo "    -o OUTPUT : output folder"
    echo "    [-h] : help"
    echo "    [-v] : version"
    exit
}

version()
{
    echo "$SOFT version $VERSION"
    exit
}

if [ $# -eq 0 ]; then
    usage
    exit
fi

while [ $# -gt 0 ]
do
    case "$1" in
	(-c) CONF=$2; shift;;
	(-i) INPUT=$2; shift;;
	(-s) STEP=$2; shift;;
	(-n) SAMPLE_NUMBER=$2; shift;;
	(-o) OUTPUT=$2; shift;;
	(-h) help;;
	(-v) version;;
	(--) shift; break;;
	(-*) die "$0: error - unrecognized option $1" 1>&2;;
	(*) break;;
    esac
    shift
done

if [[ -z ${CONF} || -z ${INPUT} || -z ${STEP} || -z ${SAMPLE_NUMBER} || -z ${OUTPUT} ]]; then
    usage
    exit
fi

##  Read configuration file  ##
source ${CONF}

##  Output folder  ##
if [ ! -d ${OUTPUT} ]; then
    mkdir -p ${OUTPUT}
fi

if [ ! -d ${OUTPUT}/${SAMPLE_NUMBER} ]; then
    mkdir -p ${OUTPUT}/${SAMPLE_NUMBER}
fi

##  STEP option  ##
for i in `echo ${STEP} | tr "," " " `
do
    NAME_STEP=${i}


##  logs folder  ##
    LOGS=${OUTPUT}/${SAMPLE_NUMBER}/logs
    if [ ! -d ${LOGS} ]; then
	mkdir -p ${LOGS}
    fi

################################
###  Clean fastq + Trimming  ###
################################

## clean_fastq step removes reads which have had a problem with the sequencer
## trimming step removes adapters in 5' and 3' side of reads
 
    OUTPUT_CLEAN=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${INPUT} | sed -e "s/.fastq.gz$/_clean.fastq/"`
    OUTPUT_TRIMMING=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_CLEAN} | sed -e 's/_clean.fastq$/_trimming_filter.fastq/'`
    if [[ ${NAME_STEP} == "trimming" || ${NAME_STEP} == "all" ]]; then
	echo "Clean fastq ..."
	cmd="${PYTHON_PATH}/python ${SCRIPTS}/selectFile.py -input ${INPUT} -output ${OUTPUT_CLEAN} 2>&1 > ${LOGS}/clean_fastq.log"
	evalecho "$cmd"

	## File verification ##
	if [ ! -f ${OUTPUT_CLEAN} ]; then
	    die " error: the *_clean.fastq file doesn't exist" 1>&2
	fi
	echo "Trimmed data ..."
	cmd="${CUTADAPT_PATH}/cutadapt -g ${ADAPTER5} -a ${ADAPTER3} -n 2 -m ${MIN_LENGTH_READS} -o ${OUTPUT_TRIMMING} ${OUTPUT_CLEAN} 2>&1 > ${LOGS}/cutadapt.log "
	evalecho "$cmd"
    fi

################
###  FASTQC  ###
################

    if [[ ${NAME_STEP} == "fastqc" || ${NAME_STEP} == "all" ]]; then
	if [ ! -f ${OUTPUT_TRIMMING} ]; then
	    die " error: the *_trimming_filter.fastq file doesn't exist" 1>&2
	fi
	echo "Fastqc ..."
        cmd="${FASTQC_PATH}/fastqc ${OUTPUT_TRIMMING}"
        evalecho "$cmd"
    fi

#################
###  Mapping  ###
#################

## In this step, we use Tophat2. We decided to apply the MAPQ.

    MAPPING=${OUTPUT}/${SAMPLE_NUMBER}/mapping_Tophat
    OUTPUT_SORT=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_TRIMMING} | sed -e 's/_trimming_filter.fastq$/_MAPQ_sort.bam/'`
    if [[ ${NAME_STEP} == "mapping" || ${NAME_STEP} == "all" ]]; then
        if [ ! -d ${MAPPING} ]; then
	    mkdir -p ${MAPPING}
    	fi

	if [ ! -f ${OUTPUT_TRIMMING} ]; then
	    die " error: the *_trimming_filter.fastq file doesn't exist" 1>&2
	fi
	echo "Mapping with Tophat2 ..."
	cmd="${TOPHAT_PATH}/tophat -o ${MAPPING} -N ${READ_MISMATCHES} --read-edit-dist ${READ_MISMATCHES} -G ${GTF_FILE} -x ${TRANSCRIPTOME_MAX_HITS} -g ${MAX_MULTIHITS} ${BOWTIE_INDEX} ${OUTPUT_TRIMMING}"
	evalecho "$cmd"

	if [ ! -f ${MAPPING}/accepted_hits.bam ]; then
	    die " error: the accepted_hits.bam file in mapping_Tophat folder doesn't exist" 1>&2
	fi
	echo "Mapping Quality >= ${MIN_MAPQ} ..."
	echo "Sorting ..."
	cmd="${SAMTOOLS_PATH}/samtools view -b -q ${MIN_MAPQ} ${MAPPING}/accepted_hits.bam | ${SAMTOOLS_PATH}/samtools sort -O bam -T ${OUTPUT}/${SAMPLE}/prefix.bam - -o ${OUTPUT_SORT}"
	evalecho "$cmd"
    fi

################################
###  Remove Duplicate reads  ###
################################

# rm_dup step permits to remove or only mark duplicated reads (Duplicated reads is a reads with the same sequence and the same length).

    OUTPUT_RM_DUPLICATED_READS=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_SORT} | sed -e 's/_MAPQ_sort.bam$/_rm_duplicated_reads.bam/'`
    OUTPUT_BAM_TO_BED=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_RM_DUPLICATED_READS} | sed -e 's/.bam$/.bed/'`

    if [[ ${NAME_STEP} == "rm_dup" || ${NAME_STEP} == "all" ]]; then
	if [ ! -f ${OUTPUT_SORT} ]; then
	    die " error: the *_MAPQ_sort.bam file doesn't exist" 1>&2
	fi
	if [ "${REMOVE_DUPLICATES}" = 1 ]; then
	    echo "Remove Duplicated reads ..."
	    echo "${OUTPUT_SORT}"
	    cmd="${JAVA_PATH}/java -Xmx30G -jar ${PICARD_TOOLS}/MarkDuplicates.jar I=${OUTPUT_SORT} O=${OUTPUT_RM_DUPLICATED_READS} METRICS_FILE=${OUTPUT}/${SAMPLE_NUMBER}/metrix.txt REMOVE_DUPLICATES=TRUE"
	    evalecho "$cmd"
	else
	    echo "Don't remove Duplicated reads"
	    cmd="${JAVA_PATH}/java -jar ${PICARD_TOOLS}/MarkDuplicates.jar I=${OUTPUT_SORT} O=${OUTPUT_RM_DUPLICATED_READS} METRICS_FILE=${OUTPUT}/${SAMPLE_NUMBER}/metrix.txt REMOVE_DUPLICATES=FALSE"
	    evalecho "$cmd"
	fi

	## Index file ##
	echo "Create a index file ..."
	cmd="${SAMTOOLS_PATH}/samtools index ${OUTPUT_RM_DUPLICATED_READS}"
	evalecho "$cmd"

       ## bamTobed ##
	if [ ! -f ${OUTPUT_RM_DUPLICATED_READS} ]; then
	    die " error: the *_rm_duplicated_reads.bam file doesn't exist" 1>&2
	fi
	cmd="${BEDTOOLS_PATH}/bamToBed -i ${OUTPUT_RM_DUPLICATED_READS} > ${OUTPUT_BAM_TO_BED}"
	evalecho "$cmd"
    fi

######################
###  Peak calling  ###
######################

## peak_calling step permits to detect peak.

    OUTPUT_PEAK_CALLING=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_BAM_TO_BED} | sed -e 's/_rm_duplicated_reads.bed$/_peaks.bed/'`

    if [[ ${NAME_STEP} == "peak_calling" || ${NAME_STEP} == "all" ]]; then
	if [ ! -f ${OUTPUT_BAM_TO_BED} ]; then
	    die " error: the *_rm_duplicated_reads.bed file doesn't exist" 1>&2
	fi
	echo "Peak calling ..."
	cmd="${PIRANHA_PATH}/Piranha ${OUTPUT_BAM_TO_BED} -s -b ${BIN_SIZE} -a ${BACKGROUND_THRESHOLD} -o ${OUTPUT_PEAK_CALLING}"
	evalecho "$cmd"
    fi

####################
###  Annotation  ###
####################

# In this step we annotate genes and after that peaks.

    OUTPUT_ADAPT_FILE=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_PEAK_CALLING} | sed -e 's/.bed$/_intersect.bed/'`
    OUTPUT_INTERSECTBED=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_ADAPT_FILE} | sed -e 's/_peaks_intersect.bed$/_intersectBED.bed/'`
    INPUT_INTERSECT_ANNOTATION=${ANNOTATION_DIR}/${ORG}/`basename whole_gene.bed`

    if [[ ${NAME_STEP} == "annotation" || ${NAME_STEP} == "all" ]]; then
	if [ ${BUILD_ANNOT} == 1 ]; then
	    echo "Build annotation file ..."             
	## Create gene file
	    cmd="${R_PATH}/R CMD BATCH \"--args infile='${ANNOTATION_DIR}/${ORG}/${UCSC_EXPORT}' out_dir='${ANNOTATION_DIR}/${ORG}/'\" ${SCRIPTS}/make_annot_gene.R ${LOGS}/make_annot_gene.Rout"
	    evalecho "$cmd"
	fi

	if [ ! -f ${OUTPUT_PEAK_CALLING} ]; then
	    die " error: the *_peaks.bed file doesn't exist" 1>&2
	fi
	echo "Treatment of peaks file to annotate peaks..."
	cmd="${R_PATH}/R CMD BATCH \"--args peakFile='${OUTPUT_PEAK_CALLING}'\" ${SCRIPTS}/adapt_file.R ${LOGS}/adapt_file.Rout"
	evalecho "$cmd"
    
	if [ ! -f ${OUTPUT_ADAPT_FILE} ]; then
	    die " error: the *_peaks_intersect.bed file doesn't exist" 1>&2
	fi
	if [ ! -f ${ANNOTATION_DIR}/${ORG}/`basename whole_gene.bed` ]; then
	    die " error: the whole_gene.bed file in annotation folder doesn't exist. Please build annotation file (BUILD_ANNOT = 1)." 1>&2
	fi
	echo "Annotation of peaks thanks to annotation of genes ..."
	cmd="${BEDTOOLS_PATH}/intersectBed -a ${OUTPUT_ADAPT_FILE} -b ${INPUT_INTERSECT_ANNOTATION} -wa -wb -s > ${OUTPUT_INTERSECTBED}"
	evalecho "$cmd"
    fi

#############################
###  Specific Annotation  ###
#############################

# In this step we annotate peaks but we choose a specific window (last exon, downstream & upstream region of peaks, peaks with downstream region).

##  Last exon  ##
    INPUT_INTERSECT_ANNOTATION_LAST_EXON=${ANNOTATION_DIR}/${ORG}/`basename last_exon_gene.bed`
    OUTPUT_INTERSECTBED_LAST_EXON=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_ADAPT_FILE} | sed -e 's/_peaks_intersect.bed$/_intersectBED_last_exon.bed/'`
    
    if [[ ${NAME_STEP} == "annot_spe" || ${NAME_STEP} == "all" ]]; then
	if [ ! -f ${ANNOTATION_DIR}/${ORG}/`basename last_exon_gene.bed` ]; then
	    die " error: the last_exon_gene.bed file in annotation folder doesn't exist. Please build annotation file (BUILD_ANNOT = 1)." 1>&2
	fi
	echo "Annotation of peaks present in the last exon region..."
	cmd="${BEDTOOLS_PATH}/intersectBed -a ${OUTPUT_ADAPT_FILE} -b ${INPUT_INTERSECT_ANNOTATION_LAST_EXON} -wa -wb -s > ${OUTPUT_INTERSECTBED_LAST_EXON}"
	evalecho "$cmd"

##  Downstream region of peak  ##
    INPUT_INTERSECT_ANNOTATION_DOWNSTREAM_END=${ANNOTATION_DIR}/${ORG}/`basename ${INPUT_INTERSECT_ANNOTATION} | sed -e 's/.bed$/_Preparative_DownstreamEnd_\${WINDOW_SIZE}_bp.bed/'`
    PRE_OUTPUT_INTERSECTBED_DOWNSTREAM=${ANNOTATION_DIR}/${ORG}/`basename ${INPUT_INTERSECT_ANNOTATION} | sed -e 's/.bed$/_\${WINDOW_SIZE}_bp_DownstreamEnd.bed/'`
    OUTPUT_INTERSECTBED_DOWNSTREAM=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_ADAPT_FILE} | sed -e 's/_peaks_intersect.bed$/_intersectBED_\${WINDOW_SIZE}_bp_DownstreamEnd.bed/'`

	if [ ${WINDOW_SIZE} != 0 ]; then
	    echo "Annotation of peaks present in the region of ${WINDOW_SIZE} bp downstream at the end of the gene ..."
	    cmd="${AWK_PATH}/awk -F\"\t\" -v WINDOW_SIZE=${WINDOW_SIZE} 'BEGIN{OFS=\"\t\"}{if(\$4~/_-/){print \$1,\$2-WINDOW_SIZE,\$2,\$4,\$5,\$6}else{print \$1,\$3,\$3+WINDOW_SIZE,\$4,\$5,\$6}}' ${INPUT_INTERSECT_ANNOTATION} > ${INPUT_INTERSECT_ANNOTATION_DOWNSTREAM_END}"
	    evalecho "$cmd"
	    cmd="${BEDTOOLS_PATH}/intersectBed -a ${INPUT_INTERSECT_ANNOTATION_DOWNSTREAM_END} -b ${INPUT_INTERSECT_ANNOTATION} -wa -s -v > ${PRE_OUTPUT_INTERSECTBED_DOWNSTREAM}"
	    evalecho "$cmd"    
	    if [ ! -f ${ANNOTATION_DIR}/${ORG}/`basename whole_gene_${WINDOW_SIZE}_bp_DownstreamEnd.bed` ]; then
		die " error: the *_${WINDOW_SIZE}_bp_DownstreamEnd.bed file doesn't exist" 1>&2
	    fi
	    cmd="${BEDTOOLS_PATH}/intersectBed -a ${OUTPUT_ADAPT_FILE} -b ${PRE_OUTPUT_INTERSECTBED_DOWNSTREAM} -wa -wb -s > ${OUTPUT_INTERSECTBED_DOWNSTREAM}"
	    evalecho "$cmd"
	fi

##  Upstream region of peak  ##
    INPUT_INTERSECT_ANNOTATION_UPSTREAM_END=${ANNOTATION_DIR}/${ORG}/`basename ${INPUT_INTERSECT_ANNOTATION} | sed -e 's/.bed$/_\${WINDOW_SIZE}_bp_UpstreamEnd.bed/'`
    OUTPUT_INTERSECTBED_UPSTREAM=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_ADAPT_FILE} | sed -e 's/_peaks_intersect.bed$/_intersectBED_\${WINDOW_SIZE}_bp_UpstreamEnd.bed/'`

	if [ ${WINDOW_SIZE} != 0 ]; then
	    echo "Annotation of peaks present in the region of ${WINDOW_SIZE} bp upstream at the end of gene ..."
	    cmd="${AWK_PATH}/awk -F\"\t\" -v WINDOW_SIZE=${WINDOW_SIZE} 'BEGIN{OFS=\"\t\"}{if(\$4~/_-/){print \$1,\$2,\$2+WINDOW_SIZE,\$4,\$5,\$6}else{print \$1,\$3-WINDOW_SIZE,\$3,\$4,\$5,\$6}}' ${INPUT_INTERSECT_ANNOTATION} > ${INPUT_INTERSECT_ANNOTATION_UPSTREAM_END}"
	    evalecho "$cmd"
	    if [ ! -f ${ANNOTATION_DIR}/${ORG}/`basename whole_gene_${WINDOW_SIZE}_bp_UpstreamEnd.bed` ]; then
		die " error: the *_${WINDOW_SIZE}_bp_UpstreamEnd.bed file doesn't exist" 1>&2
	    fi
	    cmd="${BEDTOOLS_PATH}/intersectBed -a ${OUTPUT_ADAPT_FILE} -b ${INPUT_INTERSECT_ANNOTATION_UPSTREAM_END} -wa -wb -s > ${OUTPUT_INTERSECTBED_UPSTREAM}"
	    evalecho "$cmd"
	fi


##  Peak + downstream region of peak  ##
    INPUT_INTERSECT_ANNOTATION_EXTENDED_END=${ANNOTATION_DIR}/${ORG}/`basename ${INPUT_INTERSECT_ANNOTATION} | sed -e 's/.bed$/_extended_\${WINDOW_SIZE}_bp.bed/' `
    OUTPUT_INTERSECTBED_EXTENDED=${OUTPUT}/${SAMPLE_NUMBER}/`basename ${OUTPUT_ADAPT_FILE} | sed -e 's/_peaks_intersect.bed$/_intersectBED_extended_\${WINDOW_SIZE}_bp.bed/'`


	if [ ${WINDOW_SIZE} != 0 ]; then
	    echo "Create a new annotation gene file ..."
	    cmd="${PYTHON_PATH}/python ${SCRIPTS}/extended_gene.py -inputA ${INPUT_INTERSECT_ANNOTATION} -inputB ${PRE_OUTPUT_INTERSECTBED_DOWNSTREAM} -output ${INPUT_INTERSECT_ANNOTATION_EXTENDED_END}"
	    evalecho "$cmd"
	    if [ ! -f ${ANNOTATION_DIR}/${ORG}/`basename whole_gene_extended_${WINDOW_SIZE}_bp.bed` ]; then
		die " error: the *_extended_${WINDOW_SIZE}_bp.bed file doesn't exist" 1>&2
	    fi
	    echo "Annotation of peaks present in all gene and region of ${WINDOW_SIZE} bp downstream at the end of gene ..."
	    cmd="${BEDTOOLS_PATH}/intersectBed -a ${OUTPUT_ADAPT_FILE} -b ${INPUT_INTERSECT_ANNOTATION_EXTENDED_END} -wa -wb -s > ${OUTPUT_INTERSECTBED_EXTENDED}"
	    evalecho "$cmd"
	fi
    fi

done

