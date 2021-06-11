#!/bin/bash 

function do_fastqc_trimmomatic_bismark_align_deduplicate_methylation_extraction () {
	basedir="/data/NHLBI_BCB/Thein_Lab/08_methylseq_redo"
	datadir="$basedir/01-fastqs"
	outdir="$basedir/02-methylseq-analysis-pipeline"
	bisulfite_converted_reference="/data/NHLBI_BCB/bin/hg38_Bisulfite_Genome"
	numcpus=2
	
	while read SAMPLE READ1 READ2; do
	    	echo $SAMPLE
#		mkdir -p "$outdir/$SAMPLE"
		sbatch  --job-name="${SAMPLE}" \
			--partition=norm \
			--time=5:00:00 \
			--mem=10g \
			--gres=lscratch:800 \
			--cpus-per-task=$numcpus \
			--error="$outdir/$SAMPLE.slurm.methylseq.err.txt" \
			--output="$outdir/$SAMPLE.slurm.methylseq.out.txt" \
			"$outdir/trimmomatic-fastqc-bismark-align-deduplicate-methylation-extraction.sh ${SAMPLE} $datadir $READ1 $READ2 $bisulfite_converted_reference $outdir $numcpus"
	done < "/data/NHLBI_BCB/Thein_Lab/08_methylseq_redo/sids.txt"
}

do_fastqc_trimmomatic_bismark_align_deduplicate_methylation_extraction

