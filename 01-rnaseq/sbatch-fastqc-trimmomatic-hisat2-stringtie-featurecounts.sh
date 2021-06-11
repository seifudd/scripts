#!/bin/bash 

function do_fastqc_trimmomatic_hisat_stringtie () {
	########################################################################################################################
	basedir="/data/NHLBI_BCB/Levine_Stew_Lab/06_Will_RNA-seq"
	datadir="$basedir/01-fastqs"
	outdir="$basedir/02-fastqc-trimmomatic-hisat2-featurecounts"
	reference="/data/NHLBI_BCB/bin/HISAT2-reference-genomes/grch38_tran"
	numcpus=8

	while read SAMPLE SAMPLENUM READ1 READ2; do
	    	echo ${SAMPLE}
#		mkdir "$outdir/${SAMPLE}"
		sbatch  --job-name="${SAMPLE}" \
			--partition=norm \
			--time=24:00:00 \
			--mem=64g \
			--gres=lscratch:150 \
			--cpus-per-task=$numcpus \
			--error="$outdir/${SAMPLE}.slurm.err.txt" \
			--output="$outdir/${SAMPLE}.slurm.out.txt" \
			"$outdir/fastqc-trimmomatic-hisat2-stringtie-featurecounts.sh ${SAMPLE} $datadir $READ1 $READ2 $reference $outdir $numcpus"
	done < "$basedir/sids.txt"
	########################################################################################################################
}

do_fastqc_trimmomatic_hisat_stringtie

