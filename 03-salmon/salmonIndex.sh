#!/bin/bash

mkdir -p /scratch/sack_APA/08-salmonV2
cd /scratch/sack_APA/08-salmonV2

	sbatch --cpus-per-task=8 --time=16:00:00 --mem=16g \
               -o salmonIndex.log \
	       sbatch-salmonIndex.sh 
