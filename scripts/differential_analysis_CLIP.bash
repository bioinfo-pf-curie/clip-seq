#!/bin/bash

set -o pipefail  # trace ERR through pipes
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value

SOFT="CLIP-seq"
VERSION="2.0.1"

################
## Function

function evalecho {
    echo 
#    echo $1
    echo
    eval $1
}
################

##  Initialization  ##
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
    echo "usage: $0 -c CONF -l INPUT_LIST -b BAM_RNA_FILES -s STEP -o OUTPUT"
    echo "Use -h option for more information"
    exit
}

help()
{
    echo " *** $SOFT  $VERSION ***"
    echo
    echo "OPTIONS"
    echo "    -c CONFIG : configuration file for $SOFT processing"
    echo "    -l LIST : input list file; tabulated columns with the sample id, path of bam file, condition of samples to compare and the id of condition (1 for test, 0 for ctrl)"
    echo "    -b BAM : bam RNA file; tabulated columns with the sample id (CLIP), path of RNA bam file and condition"
    echo "    -s STEP : run all or only a subset of the $SOFT workflow"
    echo " 	  merge_peaks : Creation of common peaks"
    echo "	  toc_peaks : Creation of table of counts for CLIP experiment"
    echo " 	  toc_genes : Creation of table of counts for RNA experiment"
    echo "	  cat_peaks_genes : Concatenation of tables of counts (CLIP and RNA)"
    echo " 	  differential_analysis : Make a differential analysis with samples to compare"
    echo "	  color_peaks : Creation of a file with color information. Up and Down regulated peaks (padj <0.05) are colored on red and blue respectively"
    echo "	  all : Launch all previous steps"
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
	(-l) INPUT_LIST=$2; shift;;
	(-b) BAM_RNA_FILES=$2; shift;;
	(-s) STEP=$2; shift;;
	(-o) OUTPUT=$2; shift;;
	(-h) help;;
	(-v) version;;
	(--) shift; break;;
	(-*) die "$0: error - unrecognized option $1" 1>&2;;
	(*) break;;
    esac
    shift
done

if [[ -z ${CONF} || -z ${INPUT_LIST} || -z ${STEP} || -z ${BAM_RNA_FILES} || -z ${OUTPUT} ]]; then
    usage
    exit
fi

##  Read configuration file  ##
source ${CONF}


##  Output folder  ##
if [ ! -d ${OUTPUT} ]; then
    mkdir -p ${OUTPUT}
fi


##  Output for all samples comparison  ##
odir=`echo ${COMBINE_SAMPLE} | tr "," "_"`
OUTPUT_ALL=${OUTPUT}/${odir}
mkdir -p ${OUTPUT_ALL} || die "Cannot create output folder"


## Read sample plan and put it in a hash table
declare -A SAMPLES
declare -A SAMPLES_NAME
while read line
do 
	id=`echo ${line} | ${AWK_PATH}/awk '{print $1}'`
	f=`echo ${line} | ${AWK_PATH}/awk '{print $2}'`
	n=`echo ${line} | ${AWK_PATH}/awk '{print $3}'`
	SAMPLES[${id}]="${f}"
	SAMPLES_NAME[${id}]="${n}"
done < ${INPUT_LIST}

##  STEP option  ##
for i in `echo ${STEP} | tr "," " " `
do
    NAME_STEP=${i}

### CLIP-seq
#####################
###  Merge peaks  ###
#####################

## For this step, we merge all peaks of CLIP-seq to create common peaks between all files
# bedtools version 2.21.0
    name=`echo $COMB | tr "," " "`
    if [[ ${NAME_STEP} == "merge_peaks" || ${NAME_STEP} == "all" ]]; then
	echo "Concatenate peak files ..."
	if [ -e ${OUTPUT_ALL}/all_peaks.bed ]; then
	    rm ${OUTPUT_ALL}/all_peaks.bed
	fi
	for i in `echo $COMB | tr "," " "`
	do
	    fn=${SAMPLES["${i}"]}
	    input="${fn}"
	    cat ${input} >> ${OUTPUT_ALL}/all_peaks.bed
	done

	echo "Merge Peaks ..."
	cmd="${BEDTOOLS_PATH}/sortBed -i ${OUTPUT_ALL}/all_peaks.bed | ${BEDTOOLS_PATH}/bedtools merge -s -c 10 -o collapse -delim \"|\" -i - | ${AWK_PATH}/awk  'BEGIN{OFS=\"\t\";c=1}{print \$1,\$2,\$3,\"mpeak_\"c,\$4,\$5;c=c+1}' > ${OUTPUT_ALL}/merged_peaks.bed"
	evalecho "$cmd"

       ## Keep one exemplary of gene name and remove peaks present in many genes
	cmd="${AWK_PATH}/awk '{name=\"\";p=\"T\";split(\$5,a,\"|\");for(x in a){if(a[x] ~ /^chr/){split(a[x],b,\"_\");if(name==\"\"){name=b[4];}else{if(p==\"T\" && name != b[4]){p=\"F\";}}}}; if(p==\"T\"){OFS=\"\t\"}{if(b[7]!=0){print \$1,\$2,\$3,\$3-\$2+1,\$3-\$2+1,b[7],\$4\"|\"\$1\"_\"b[2]\"_\"b[3]\"_\"name\"_\"b[5]\"_\"b[6]\"_\"b[7],\$1\"_\"b[2]\"_\"b[3]\"_\"name\"_\"b[6];}else{print \$1,\$2,\$3,\$3-\$2+1,\$3-\$2+1,b[6],\$4\"|\"\$1\"_\"b[2]\"_\"b[3]\"_\"name\"_\"b[5]\"_\"b[6],\$1\"_\"b[2]\"_\"b[3]\"_\"name;}}}' ${OUTPUT_ALL}/merged_peaks.bed > ${OUTPUT_ALL}/merged_peaks_finallist.bed"
	evalecho "$cmd"

       ## Add information: CLIP| for NR peak and CLIP for NM
	cmd="${AWK_PATH}/awk 'BEGIN{OFS=\"\t\"}{if(\$7 ~/_NR_/){print \$0,\"CLIP|\"}else{print \$0,\"CLIP\"}}' ${OUTPUT_ALL}/merged_peaks_finallist.bed > ${OUTPUT_ALL}/merged_peaks_finallist_annot.bed"
	evalecho "$cmd"
    fi

#########################################
###  Creation of the table of counts  ###
#########################################
    name=`echo $COMB | tr "," "\t"`
    if [[ ${NAME_STEP} == "toc_peaks" || ${NAME_STEP} == "all" ]]; then
	echo "Creation table of counts with CLIP-seq data..."
	for i in `echo ${COMB} | tr "," " "`
	do
	    fn=${SAMPLES["${i}"]}
	    dname=`basename ${fn} | sed -e 's/_intersectBED.*/_MAPQ_rm_duplicated_reads.bam/'`
	    file="${file} ${OUTPUT}/${i}/${dname}"
	done
	cmd="${BEDTOOLS_PATH}/bedtools multicov -s -bams ${file} -bed ${OUTPUT_ALL}/merged_peaks_finallist_annot.bed | ${AWK_PATH}/awk -v name=\"\$name\" 'BEGIN{OFS=\"\t\";print \"chr\", \"start\", \"end\", \"length\", \"length2\", \"strand\", \"peak\", \"gene\", \"status\", name}{print \$0}' > ${OUTPUT_ALL}/merged_peaks_finallist_annot_qt.bed"
	evalecho "$cmd"

       ## File improvement
	cmd="${PYTHON_PATH}/python ${SCRIPTS}/sort_peak.py -i ${OUTPUT_ALL}/merged_peaks_finallist_annot_qt.bed -o ${OUTPUT_ALL}/merged_peaks_finallist_IanDif.bed"
	evalecho "$cmd"
    fi


### RNA-seq
#####################################################
###  Read sample plan and put it in a hash table  ###
#####################################################
    declare -A SAMPLES_RNA
    declare -A SAMPLES_RNA_NAME
    while read line
    do 
	idR=`echo ${line} | ${AWK_PATH}/awk '{print $1}'`
	fR=`echo ${line} | ${AWK_PATH}/awk '{print $2}'`
	nR=`echo ${line} | ${AWK_PATH}/awk '{print $3}'`
	SAMPLES_RNA[${idR}]="${fR}"
	SAMPLES_RNA_NAME[${idR}]="${nR}"
	fileR="${fileR} ${fR}"
    done < ${BAM_RNA_FILES}

#######################################################
###  Create table of counts with R (featureCounts)  ###
#######################################################

    if [[ ${NAME_STEP} == "toc_genes" || ${NAME_STEP} == "all" ]]; then
	echo "Creation table of counts with RNA-seq data..."
	cmd="${R_PATH}/Rscript "${SCRIPTS}"/count_RNA-seq.R --args "${ANNOTATION_DIR}"/"${ORG}"/whole_gene.bed \""${fileR}'" '${OUTPUT_ALL}"/countTable_RNA.bed > "${LOGS}"/count_RNA-seq.Rout"
	evalecho "$cmd"

       #Add information RNA| if NR gene and RNA else
	cmd="${AWK_PATH}/awk 'BEGIN{OFS=\"\t\"}{if(\$1 ~/_NR_/){print \$0,\"RNA|\"}else{print \$0,\"RNA\"}}' ${OUTPUT_ALL}/countTable_RNA.bed > ${OUTPUT_ALL}/whole_gene_finallist.bed"
	evalecho "$cmd"

       #Add geneName column (chr_start_end_gene)
	cmd="${PYTHON_PATH}/python ${SCRIPTS}/sort_gene.py -i ${OUTPUT_ALL}/whole_gene_finallist.bed -o ${OUTPUT_ALL}/whole_gene_finallist_IanDif.bed"
	evalecho "$cmd"
    fi

#######################################
##  Cat CLIP-seq and RNA-seq files  ###
#######################################

    if [[ ${NAME_STEP} == "cat_peaks_genes" || ${NAME_STEP} == "all" ]]; then
	echo "Concatenate CLIP-seq file and RNA-seq file ..."
	if [ -e ${OUTPUT_ALL}/all_files.bed ]; then
	    rm ${OUTPUT_ALL}/all_files.bed
	fi
	for i in `echo $COMB | tr "," " "`
	do
	    cat ${OUTPUT_ALL}/merged_peaks_finallist_IanDif.bed ${OUTPUT_ALL}/whole_gene_finallist_IanDif.bed > ${OUTPUT_ALL}/all_files.bed
	done
    fi

###########################################
###  Differential analysis with DESeq2  ###
###########################################
    if [[ ${NAME_STEP} == "differential_analysis" || ${NAME_STEP} == "all" ]]; then
	echo "Differential analysis between CLIP-seq and RNA-seq ..."
	cmd="${R_PATH}/R CMD BATCH \"--args peakfile='${OUTPUT_ALL}/all_files.bed' input_list='${INPUT_LIST}' polyA_lib='${SCRIPTS}/polyA_lib.R' min_count_cond='${MIN_COUNT_PER_COND}'\" ${SCRIPTS}/compare.R ${LOGS}/compare.Rout"
	evalecho "$cmd"
    fi

###################
###  Add color  ###
###################

## Add color information. up & down regulated peaks (determined thanks to padj & log2FC_CLIP-log2FC_RNA) are in red and blue respectively.

	if [[ ${NAME_STEP} == "color_peaks" || ${NAME_STEP} == "all" ]]; then
		echo "Add color for differential peaks ..."
		cmd="${PYTHON_PATH}/python ${SCRIPTS}/add_color_padj.py -input ${OUTPUT_ALL}/all_files_peakClip_SumRNA_total.txt -output ${OUTPUT_ALL}/all_peakCLIP_color.txt"
		evalecho "$cmd"
	fi
	
done
