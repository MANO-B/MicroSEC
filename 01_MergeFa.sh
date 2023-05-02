#!/bin/bash

# Filename : $RCSfile: 01_MergeFa.sh,v $
# Version  : $Revision: 1.9 $
# Updated  : $Date: 2020-10-19 16:59:06+09 $
# Copyright(C) National Cancer Center Japan 2014

# Unzip and concatenate fasta gzip file(s) into 1 file.

set -ue;

EXE_D=$( cd $( dirname $0 ); pwd )
DATA_D=$( dirname "${EXE_D}" )/data
IN_D=${DATA_D}/00_org

OUT_D="${DATA_D}"/$( basename -s .sh $0 )
mkdir -p "${OUT_D}";

OUT_F="${OUT_D}/rna.fa"
LOG_F="${OUT_D}"/log_$( basename -s .sh $0 ).txt

echo -n > ${OUT_F};
echo -n > ${LOG_F};

for IN_F in ${IN_D}/*.gz
do
    echo ${IN_F};
    gzip -dc "${IN_F}" >> "${OUT_F}";

    echo "${IN_F}" >> "${LOG_F}";
    zgrep "^>" "${IN_F}" | sed -e 's/_.*$//' | sort | uniq -c >> "${LOG_F}";
done
