#!/bin/bash

set -e

################################################################################
# Set Variables                                                                #
################################################################################

trimm_th=10
trimm_rm=1
#trimm_rm=0
#index_builder=1
#db_ftp="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/latest_release/"
#transcriptome="gencode.*.transcripts.fa.gz"
#assembly="GRCh38.primary_assembly.genome.fa.gz"
index_rm=0
index_builder=0
index_path="./Index/index"
salmon_th=20

################################################################################
################################################################################
# Main Program                                                                 #
################################################################################
################################################################################

source activate fastq_to_transcripts

################################################################################
# Salmon Index                                                                 #
################################################################################

if [ $index_builder == 1 ]; then
  echo "Index Builder Mode Enabled: Creating New Index"

  mkdir Index

  wget -r -l1 -np -nd -A $transcriptome -P ./Index/ $db_ftp
  wget -r -l1 -np -nd -A $assembly -P ./Index/ $db_ftp

  grep "^>" <(gunzip -c ./Index/$assembly) | cut -d " " -f 1 > \
    ./Index/decoys.txt
  sed -i -e 's/>//g' ./Index/decoys.txt
  cat ./Index/$transcriptome ./Index/$assembly > ./Index/gentrome.fa.gz

  salmon index \
              -t ./Index/gentrome.fa.gz -d ./Index/decoys.txt \
              -k 31 -p $salmon_th -i ./Index/index --gencode

  index_path="./Index/index"

elif [ $index_builder == 0 ]; then
  echo "Index Builder Mode Disabled: Using the Index located in $index_path"

else
  echo "Error: Failed to detect Index builder status. Please set the variable" \
    "index_builder=1 to create a new index or index_builder=0 to use one" \
    "already created"

  ls ./Index 2> /dev/null

fi

################################################################################
# QC and Trimming                                                              #
################################################################################

mkdir QC trimmed

for f1 in ./fastq/*_1.fastq.gz; do
  f2=${f1%%_1.fastq.gz}"_2.fastq.gz"
  samp1=`basename ${f1%%.fastq.gz}`
  samp2=`basename ${f2%%.fastq.gz}`
  echo "Performing QC and trimming of samples ${samp1} and ${samp2}"
  fastp -i $f1 -I $f2 \
        -o ./trimmed/${samp1}_trimmed.fastq.gz \
        -O ./trimmed/${samp2}_trimmed.fastq.gz \
        -w $trimm_th \
        --detect_adapter_for_pe --trim_poly_x --correction -r -M 10 -l 20 \
        -h ./QC/${samp1%%_1}_report.html -j ./QC/${samp1%%_1}_report.json \
        -R "Fastp Report of ${samp1%%_1} Sample"
done

################################################################################
# Transcript Quantification                                                    #
################################################################################

mkdir quants

for t1 in ./trimmed/*_1_trimmed.fastq.gz
do
  t2=${t1%%_1_trimmed.fastq.gz}"_2_trimmed.fastq.gz"
  samp=`basename ${t1%%_1_trimmed.fastq.gz}`
  echo "Performing Salmon Quantification (Mapping-based mode) of Sample ${samp}"
  salmon quant \
              -i $index_path -l A -1 ${t1} -2 ${t2} -o quants/${samp} \
              -p $salmon_th --validateMappings --seqBias --gcBias #\
              #--numBootstraps 100
done

################################################################################
# Cleaning                                                                     #
################################################################################

conda deactivate

if [ $trimm_rm == 1 ]; then
  echo "Removing trimmed fastq files"
  rm -r ./trimmed
fi

if [ $index_rm == 1 ]; then
  echo "Removing index files"
  rm -r ./Index
fi
################################################################################
