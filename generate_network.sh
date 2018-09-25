#! /bin/bash
SAMTOOLS=samtools
SAMBAMBA=sambamba
R_SCRIPT=Rscript
PYTHON=python
BEDTOOLS=bedtools
BLASTN=blastn
BLAST_MKDB=makeblastdb

SCRIPT_NAME=$(basename $0)
DIR_NAME="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "$(date +"%T")" : script begins.
IFS=$''
while getopts b:g:o:c:p:l:d:n:u:r:a:s:m:h ARG
do
  case ${ARG} in
    (b) BAM_FILE=$(readlink -f "$OPTARG");;
    (g) GTF_FILE=$(readlink -f "$OPTARG");;
    (o) OUTPUT_DIRECTORY=$(readlink -f "$OPTARG");;
    (c) CACHE_DIRECTORY="$OPTARG";;
    (p) PERCENTAGE="$OPTARG";;
    (l) COVERAGE="$OPTARG";;
    (d) GENE_LIST="$OPTARG";;
    (n) GENE_LIST=($(cat "$OPTARG"));;
    (u) UNIFICATION="$OPTARG";;
    (r) COORD_RANGE_INT="$OPTARG";;
    (a) SAMP_ANNOTATION="$OPTARG";;
    (s) SORTED_BAM_FOLDER="$OPTARG";;
    (f) FASTQ_FILE=$(readlink -f "$OPTARG");;
    (h)
      echo " -b <directory> The input folder with all BAM files"
      echo " -g <file> The input GTF file"
      echo " -o <The directory for writing the output files> directory"
      echo " -c <The cache directory> directory"
      echo " -p <The percentage similarity value (default 85)> value"
      echo " -l <The percentage coverage value (default 30)> value"
      echo " -d <A list of genes to examine ("GENE1, GENE2, ...")> list"
      echo " -n <file_array> A text file with gene names in each line"
      echo " -u <The method used to unify reads. Options: by sequence (by_seq, default), by coordinates (by_coord) or by binned chromosome ranges (no_blast) > value"
      echo " -r <The interval lengths used in unification using coordinate range> value"
      echo " -a <If the file is supplied, it is used to average reads within specified groups> value"
      echo " -s <If this value is supplied, skip step that sorts/indexes .bam files> value"
      exit 0;;
    (*)
      exit 1;;
  esac
done

if [ ! "${BAM_FILE}" ] || [ ! -e "${BAM_FILE}" ];
then
	echo "BAM file not supplied or not found; use -b option"
	exit 1;
fi

if [ ! "${GTF_FILE}" ] || [ ! -e "${GTF_FILE}" ];
then
	echo "GTF file not supplied or not found; use -g option"
	exit 1;
fi

if [ ! "${OUTPUT_DIRECTORY}" ];
then
	OUTPUT_DIRECTORY=$(readlink -f output)
fi

if [ ! "${CACHE_DIRECTORY}" ];
then
	CACHE_DIRECTORY=$(readlink -f cache)
	INPUT_HASH=$(cat ${UNSORTED_BAM_FILE} ${CHROMOSOME_LENGTH_FILE} ${GTF_FILE} | md5sum)
	RANDOM_NUM1=$(($RANDOM%999))
	RANDOM_NUM2=$(($RANDOM%999))
	RANDOM_NUM3=$(($RANDOM%999))
	INPUT_HASH=${INPUT_HASH%% *}
	mkdir -p $OUTPUT_DIRECTORY
	INPUT_HASH="${INPUT_HASH}_${RANDOM_NUM1}_${RANDOM_NUM2}_${RANDOM_NUM3}"
	if [ -d ${CACHE_DIRECTORY}/${INPUT_HASH} ]; then
		RANDOM_NUM4=$(($RANDOM%999))
		INPUT_HASH="${INPUT_HASH}_${RANDOM_NUM4}"
	fi
	CACHE_DIRECTORY="${CACHE_DIRECTORY}/${INPUT_HASH}"
else
	CACHE_DIRECTORY=$(readlink -f cache)/$CACHE_DIRECTORY
fi

if [ ! "${PERCENTAGE}" ];
then
	PERCENTAGE="85"
fi

if [ ! "${COVERAGE}" ];
then
	COVERAGE="30"
fi

if [ ! "${GENE_LIST}" ];
then
	echo "Gene list not supplied or not found; use -d option"
	exit 1;
fi

if [ ! "${GENE_LIST}" ];
then
	echo "Gene list not supplied or not found; use -n option"
	exit 1;
fi

if [ ! "${UNIFICATION}" ];
then
	UNIFICATION="by_coord"
else
	UNIFICATION=$(echo "$UNIFICATION" | tr '[:upper:]' '[:lower:]')
fi

if [ ! "${COORD_RANGE_INT}" ];
then
	COORD_RANGE_INT="10"
fi

if [ ! "${SAMP_ANNOTATION}" ];
then
	SAMP_ANNOTATION=" "
fi

if [ ! "${SORTED_BAM_FOLDER}" ];
then
	SORTED_BAM_FOLDER=" "
fi


mkdir -p $CACHE_DIRECTORY

BASENAME_BAM_FILE=$(basename ${BAM_FILE})
NO_EXT_BAM_FILE="${BASENAME_BAM_FILE%.*}"

# Print out the sample list from sample annotation, if there is input sample annotation
if [ "$SAMP_ANNOTATION" != " " ];
then
	awk '{print $1}' $SAMP_ANNOTATION > $CACHE_DIRECTORY/current_samples.txt
else
	rm -f $CACHE_DIRECTORY/current_samples.txt
	for file in ${BAM_FILE}/*.bam
	do 
		tmp=${file##*/}
		echo ${tmp%.bam} >> $CACHE_DIRECTORY/current_samples.txt
	done
fi

if [ "$SORTED_BAM_FOLDER" == " " ];
then
	SORTED_BAM_FOLDER="${CACHE_DIRECTORY}/sorted"
	SORT_TMP_DIR="${CACHE_DIRECTORY}/sambamba_tmp"
	rm -rf $SORTED_BAM_FOLDER
	mkdir -p $SORTED_BAM_FOLDER
	#for FILE in $BAM_FILE/*.bam; do
	while read -r sample
	do
		FILE="${BAM_FILE}/${sample}.bam"
		BASENAME_BAM_FILE=$(basename ${FILE})
		NO_EXT_BAM_FILE="${BASENAME_BAM_FILE%.*}"
		SORTED_BAM="${SORTED_BAM_FOLDER}/${NO_EXT_BAM_FILE}.bam"
		if ${SAMBAMBA} sort ${FILE} -o ${SORTED_BAM} --tmpdir=$SORT_TMP_DIR
		then
			rm ${FILE}
			rm -r $SORT_TMP_DIR
			mv ${SORTED_BAM} ${FILE}
		fi
		${SAMBAMBA} index ${FILE}
	done < ${CACHE_DIRECTORY}/current_samples.txt
fi

VALID_GENES=$(perl -pe 's/.*gene_name "([^"]+)".*/\1/' ${GTF_FILE} | tr '[:lower:]' '[:upper:]' | sort | uniq)
GENE_LIST=$(echo "$GENE_LIST" | tr '[:lower:]' '[:upper:]' | perl -pe 's/\s*,\s*|\s+/ /g')
GENE_LIST=$(echo ${GENE_LIST// /|})
python py_identify_region.py ${GTF_FILE} ${CACHE_DIRECTORY} ${GENE_LIST}
CHROMOSOME_FILENAME="${CACHE_DIRECTORY}/chromosome_region.txt"

while read -r line
do
	CURRENT_GENE="$(echo $line | cut -d ',' -f 1)"
	CURRENT_CHR_LOC="$(echo $line | cut -d ',' -f 2)"
	CURRENT_STRAND="$(echo $line | cut -d ',' -f 3)"
	echo "$(date +"%T")" : current gene/chromosome loci/strand = $CURRENT_GENE/$CURRENT_CHR_LOC/$CURRENT_STRAND.
	#for FILE in $BAM_FILE/*.bam; do
	while read -r sample
	do
		FILE="${BAM_FILE}/${sample}.bam"
		BASENAME_BAM_FILE=$(basename ${FILE})
		NO_EXT_BAM_FILE="${BASENAME_BAM_FILE%.*}"
		FILTERED="${CACHE_DIRECTORY}/$CURRENT_GENE/"
		mkdir -p "${FILTERED}bam"
		mkdir -p "${FILTERED}sam"
		mkdir -p "${FILTERED}bed"
		mkdir -p "${FILTERED}node"
		mkdir -p "${FILTERED}fsa"
		FILTERED_BAM="${FILTERED}bam/$NO_EXT_BAM_FILE.bam"
		FILTERED_SAM="${FILTERED}sam/$NO_EXT_BAM_FILE.sam"
		FILTERED_BED="${FILTERED}bed/$NO_EXT_BAM_FILE.bed"
		${SAMBAMBA} slice ${FILE} $CURRENT_CHR_LOC -o ${FILTERED_BAM}
		${SAMBAMBA} view -h ${FILTERED_BAM} > ${FILTERED_SAM}
		${BEDTOOLS} bamtobed -i ${FILTERED_BAM} -bed12 > ${FILTERED_BED}
	done < ${CACHE_DIRECTORY}/current_samples.txt
	NUM_PAIR2=$(${SAMTOOLS} view -c -f 1 ${FILTERED_BAM})
	echo "$(date +"%T")" : bam converted to sam/bed.
	${PYTHON} py_sam2fasta.py ${CACHE_DIRECTORY} ${OUTPUT_DIRECTORY} ${UNIFICATION} ${CURRENT_GENE} ${CURRENT_STRAND} ${CURRENT_CHR_LOC} ${SAMP_ANNOTATION} ${NUM_PAIR2}
	echo "$(date +"%T")" : sam2fasta complete.
	if [[ $UNIFICATION != *"no_blast"* ]]
	then
		FSA_DB="${CACHE_DIRECTORY}/$CURRENT_GENE/fsa/$CURRENT_GENE.fsa"
		BLASTN_OUTPUT="${CACHE_DIRECTORY}/$CURRENT_GENE/blast"
		$BLAST_MKDB -in $FSA_DB -parse_seqids -dbtype nucl
		mkdir -p $BLASTN_OUTPUT
		$BLASTN -query $FSA_DB \
		-db $FSA_DB \
		-task blastn \
		-out "$BLASTN_OUTPUT/$CURRENT_GENE-blast.tab" \
		-num_threads 15 -outfmt 6 \
		-max_target_seqs 100 \
		-evalue 1e-5 \
		-perc_identity $PERCENTAGE
		${PYTHON} py_blast2layout.py ${CACHE_DIRECTORY} ${CURRENT_GENE} ${OUTPUT_DIRECTORY} $PERCENTAGE $COVERAGE ${UNIFICATION} ${COORD_RANGE_INT}
	#else
	#	Rscript r_stats_differential_alternative_splicing.r ${CACHE_DIRECTORY}/${CURRENT_GENE} "${OUTPUT_DIRECTORY}/${INPUT_HASH}/${CURRENT_GENE}_${UNIFICATION}.gml"
	fi
	#rm -rf "${CACHE_DIRECTORY}/${CURRENT_GENE}"
done < "$CHROMOSOME_FILENAME"

echo "$(date +"%T")" : script ends.
