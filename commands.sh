#!/bin/bash

j=6
for i in 2 3 4
do
	cp /data/NHLBI_BCB/NIDCR/01_healthy_scRNA/resequenced/Tm242/Tm242_S1_L00${i}_I1_001.fastq.gz /data/NHLBI_BCB/NIDCR/01_healthy_scRNA/Tm242/Tm242_S1_L00${j}_I1_001.fastq.gz
	cp /data/NHLBI_BCB/NIDCR/01_healthy_scRNA/resequenced/Tm242/Tm242_S1_L00${i}_R1_001.fastq.gz /data/NHLBI_BCB/NIDCR/01_healthy_scRNA/Tm242/Tm242_S1_L00${j}_R1_001.fastq.gz
	cp /data/NHLBI_BCB/NIDCR/01_healthy_scRNA/resequenced/Tm242/Tm242_S1_L00${i}_R2_001.fastq.gz /data/NHLBI_BCB/NIDCR/01_healthy_scRNA/Tm242/Tm242_S1_L00${j}_R2_001.fastq.gz
	let "j=j+1"
	echo $i,$j
done
