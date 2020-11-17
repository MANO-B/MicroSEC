#!/bin/sh
#####
LOG_E=$1
LOG_O=$2
SRC=$3
BIN=$4
UTIL=$5
RESULT_DIR=$6
SLOT=$7
SAMPLE=$8
CHECK_TYPE=$9
shift 9
ADAPTER_SEQ_1=$1
ADAPTER_SEQ_2=$2
GENOME=$3
#####
source $UTIL
#####
BAM=090_merge_bam
MUT_ID=020_mut_call
SUM=051_read_summary
PRE=261_mutation_list_with_peptide
CUR=265_MicroSEC
EXE="Rscript ${BIN}/step_${CUR}.R
#####
BAM_FILE=${RESULT_DIR}/${BAM}_${SAMPLE}/${SAMPLE}.realigned.bam
MUT_ID_DIR=${RESULT_DIR}/${MUT_ID}_${SAMPLE}
SUM_FILE=${RESULT_DIR}/${SUM}_${SAMPLE}/${SAMPLE}_RS.txt 
PRE_DIR=${RESULT_DIR}/${PRE}_${SAMPLE}
PRE_FINISH_FILE_NO=`ls ${PRE_DIR}/*_finish.txt | wc -l`
if [ $PRE_FINISH_FILE_NO -ne $CHECK_NO ]; then
        exit
fi

OUT_DIR=${RESULT_DIR}/${CUR}_${SAMPLE}

CUR_FINISH_FILE=${OUT_DIR}/${SAMPLE}-${CHECK_TYPE}_finish.txt
if [ -e $CUR_FINISH_FILE ]; then
        exit
fi

OUT=${OUT_DIR}/${SAMPLE}_${CHECK_TYPE}.gz
if [ ! -e $CUR_FINISH_FILE ]; then
        rm -rf $OUT
fi

LOGS=($LOG_E $LOG_O)
for log in ${LOGS[@]}
        do
        if [ -s $log ]; then
                : > $log
        fi
done

if [ ! -e $OUT_DIR ]; then
        mkdir $OUT_DIR
fi

FILE=${PRE_DIR}/${SAMPLE}_${CHECK_TYPE}.gz

READ=`awk 'NR % 2 == 0 {print $3}' $SUM_FILE`

$EXE $OUT $SAMPLE $FILE $BAM_FILE $MUT_ID_DIR $READ $ADAPTER_SEQ_1 $ADAPTER_SEQ_2 $GENOME
check_error $?

touch $CUR_FINISH_FILE
