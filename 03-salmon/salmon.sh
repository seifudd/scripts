#!/bin/bash 


mkdir -p /scratch/sack_APA/08-salmonV2
cd /scratch/sack_APA/08-salmonV2

for i in {3..63}; do
#for i in  1 2 ; do
    SAMPLE="S$i"
	echo $SAMPLE
	sbatch --cpus-per-task=8 --time=96:00:00 --mem=16g -J $SAMPLE \
               -o $SAMPLE.log \
	       10-sbatch-salmon.sh $SAMPLE  /scratch/sack_APA/01-fastqs
done

