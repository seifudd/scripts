#!/bin/bash

module load salmon


SAMPLE=$1
DATAPATH=$2


salmon quant -i genomeIndexHS -l IU -1 $DATAPATH/${SAMPLE}_R1.fastq.gz -2 $DATAPATH/${SAMPLE}_R2.fastq.gz -o ${SAMPLE}.quant 
