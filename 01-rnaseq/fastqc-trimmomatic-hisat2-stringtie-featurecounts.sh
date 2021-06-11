#!/bin/bash

set -e        # stop the script if a command fails

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

SAMPLE=$1
DATAPATH=$2
READ1=$3
READ2=$4
REF=$5
outdir=$6
numcpus=$7

module load hisat
module load samtools
module load trimmomatic
module load fastqc
module load stringtie
module load subread

function do_fastqc () {
	date
	########################################################################################################################
	# if trimming, change DATAPATH
	DATAPATH="/lscratch/${SLURM_JOBID}"

	# if trimming, change $READ1 and $READ2
	READ1="${SAMPLE}_1P.fastq.gz"
	READ2="${SAMPLE}_2P.fastq.gz"

	# if trimming, change $out_dir to something like "fastqc_post_trimming" if you prefer
#	out_dir="fastqc"
	out_dir="fastqc_post_trimming"

#	fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN

	mkdir -p $outdir/$SAMPLE/$out_dir

	fastqc -o "$outdir/$SAMPLE/$out_dir"  \
	--nogroup \
	"$DATAPATH/$READ1"  \
	"$DATAPATH/$READ2"  \
	|| fail "fastqc failed"

	echo "fastqc done"
	########################################################################################################################
	date
}

function do_trimmomatic () {
	date
	########################################################################################################################
	java -jar $TRIMMOJAR PE \
	            -threads $numcpus \
	            "$DATAPATH/$READ1" \
	            "$DATAPATH/$READ2" \
	 	    -baseout "/lscratch/${SLURM_JOBID}/${SAMPLE}.fastq.gz" \
	           ILLUMINACLIP:"/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa":2:30:10 \
	           MINLEN:50

	echo "trimmomatic done"
	########################################################################################################################
	date
}

function do_hisat2 () {
	date
	########################################################################################################################
	# if trimming, please change DATAPATH
	DATAPATH="/lscratch/${SLURM_JOBID}"

	# if trimming, please change $READ1 and $READ2
	READ1="${SAMPLE}_1P.fastq.gz"
	READ2="${SAMPLE}_2P.fastq.gz"

	out_dir="hisat2"
	mkdir -p $outdir/$SAMPLE/$out_dir

	hisat2 	-p $numcpus \
			-x $REF/genome_tran \
			--downstream-transcriptome-assembly \
			-1 "$DATAPATH/$READ1" \
			-2 "$DATAPATH/$READ2" \
			--rg-id $SAMPLE --rg SM:$SAMPLE \
	| samtools view -h -f 3 -O SAM - \
	| perl -nle  'print if m/^@(?:[A-Z]{2})\s|\bNH:i:1\b/' \
	| samtools sort -@ $numcpus \
		-o "$outdir/$SAMPLE/$out_dir/$SAMPLE.unique.bam" \
		-T /lscratch/${SLURM_JOB_ID}/${SAMPLE}_chunk -
	   
	samtools index "$outdir/$SAMPLE/$out_dir/$SAMPLE.unique.bam"
	echo "hisat2 done"
	########################################################################################################################
	date
}

function do_stringtie () {
	date
	########################################################################################################################
	GTF="/data/NHLBI_BCB/bin/HISAT2-reference-genomes/GENCODE_human_v28/gencode.v28.annotation.gtf"
	bamfile="$outdir/$SAMPLE/hisat2/$SAMPLE.unique.bam"
	out_dir="stringtie"
	mkdir -p $outdir/$SAMPLE/$out_dir

	# only reference
	stringtie -p 4 \
	  -o $outdir/$SAMPLE/$out_dir/$SAMPLE.gtf \
	  -e -G $GTF \
	  -l $SAMPLE \
	  -v -B -c 10 -j 5  -f 0.1 \
	  $bamfile

	# novel
	# stringtie -p 4 \
	#  -o $outdir/$SAMPLE/$out_dir/$SAMPLE.gtf \
	#  -G $GTF \
	#  -l $SAMPLE \
	#  -v -B -c 10 -j 5  -f 0.1 \
	#  $bamfile

	# https://github.com/gpertea/stringtie
	echo -e "
	--version : print current version at stdout
	 -h print this usage message
	 -G reference annotation to use for guiding the assembly process (GTF/GFF3)
	 -l name prefix for output transcripts (default: STRG)
	 -f minimum isoform fraction (default: 0.1)
	 -m minimum assembled transcript length to report (default 100bp)
	 -o output path/file name for the assembled transcripts GTF (default: stdout)
	 -a minimum anchor length for junctions (default: 10)
	 -j minimum junction coverage (default: 1)
	 -t disable trimming of predicted transcripts based on coverage
	    (default: coverage trimming is enabled)
	 -c minimum reads per bp coverage to consider for transcript assembly (default: 2.5)
	 -v verbose (log bundle processing details)
	 -g gap between read mappings triggering a new bundle (default: 50)
	 -C output file with reference transcripts that are covered by reads
	 -M fraction of bundle allowed to be covered by multi-hit reads (default:0.95)
	 -p number of threads (CPUs) to use (default: 1)
	 -A gene abundance estimation output file name
	 -B enable output of Ballgown table files which will be created in the
	    same directory as the output GTF (requires -G, -o recommended)
	 -b enable output of Ballgown table files but these files will be 
	    created under the directory path given as <dir_path>
	 -e only estimates the abundance of given reference transcripts (requires -G)
	 -x do not assemble any transcripts on the given reference sequence(s)

	Transcript merge usage mode:

	 stringtie --merge [Options] { gtf_list | strg1.gtf ...}
	With this option StringTie will assemble transcripts from multiple
	input files generating a unified non-redundant set of isoforms. In this
	usage mode the following options are available:
	  -G <guide_gff>   reference annotation to include in the merging (GTF/GFF3)
	  -o <out_gtf>     output file name for the merged transcripts GTF
		            (default: stdout)
	  -m <min_len>     minimum input transcript length to include in the merge
		            (default: 50)
	  -c <min_cov>     minimum input transcript coverage to include in the merge
		            (default: 0)
	  -F <min_fpkm>    minimum input transcript FPKM to include in the merge
		            (default: 1.0)
	  -T <min_tpm>     minimum input transcript TPM to include in the merge
		            (default: 1.0)
	  -f <min_iso>     minimum isoform fraction (default: 0.01)
	  -g <gap_len>     gap between transcripts to merge together (default: 250)
	  -i               keep merged transcripts with retained introns; by default
		           these are not kept unless there is strong evidence for them
	  -l <label>       name prefix for output transcripts (default: MSTRG)
	" > /dev/null
	echo "stringtie done"
	########################################################################################################################
	date
}

function do_featurecounts () {
	date
	########################################################################################################################
	GTF="/data/NHLBI_BCB/bin/HISAT2-reference-genomes/GENCODE_human_v28/gencode.v28.annotation.gtf"
	bamfile="$outdir/$SAMPLE/hisat2/$SAMPLE.unique.bam"

	out_dir="featurecounts"
	mkdir -p $outdir/$SAMPLE/$out_dir

	s=0  # -s strand-specific : 0 (unstranded), 1 (stranded) 2 (reversely stranded). 0 default.
	# -f -Q -M -O
	
	end_num="" # if single end sequencing, please change to "single"

	if [[ $end_num == "single" ]]
	then

	featureCounts -T $numcpus \
		-t exon \
		-g gene_id \
		-a $GTF \
		-s $s \
		-o $outdir/$SAMPLE/$out_dir/$SAMPLE.genefeatureCounts.txt \
		$bamfile

	else

	featureCounts -T $numcpus \
		-t exon \
		-g gene_id \
		-a $GTF \
		-s $s -p -M -O \
		-o $outdir/$SAMPLE/$out_dir/$SAMPLE.genefeatureCounts.txt \
		$bamfile

	fi

	echo -e "
	Version 1.4.6-p3

	Usage: featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ... 

	    Required parameters:

	    -a <input>	Give the name of the annotation file. The program assumes
		      	that the provided annotation file is in GTF format. Use -F
		      	option to specify other annotation formats.
	    
	    -o <input>	Give the name of the output file. The output file contains
		      	the number of reads assigned to each meta-feature (or each
		      	feature if -f is specified). A meta-feature is the aggregation
		      	of features, grouped by using gene identifiers. Please refer
		      	to the users guide for more details.
	    
	   input_files	Give the names of input read files that include the read
		      	mapping results. Format of input files is automatically
		      	determined (SAM or BAM). Paired-end reads will be
		      	automatically re-ordered if it is found that reads from the
		      	same pair are not adjacent to each other. Multiple files can
		      	be provided at the same time.
	    
	    Optional parameters:
	    
	    -A <input>	Specify the name of a file including aliases of chromosome
		      	names. The file should be a comma delimited text file that
		      	includes two columns. The first column gives the chromosome
		      	names used in the annotation and the second column gives the
		      	chromosome names used by reads. This file should not contain
		      	header lines. Names included in this file are case sensitive.
	    
	    -F <input>	Specify the format of the annotation file. Acceptable formats
		      	include 'GTF' and 'SAF'. 'GTF' by default. Please refer to the
		      	users guide for SAF annotation format.
	    
	    -t <input>	Specify the feature type. Only rows which have the matched
		      	matched feature type in the provided GTF annotation file
		      	will be included for read counting. 'exon' by default.
	    
	    -g <input>	Specify the attribute type used to group features (eg. exons)
		      	into meta-features (eg. genes), when GTF annotation is provided.
		      	'gene_id' by default. This attribute type is usually the gene
		      	identifier. This argument is useful for the meta-feature level
		      	summarization.
	    
	    -f        	If specified, read summarization will be performed at the 
		      	feature level (eg. exon level). Otherwise, it is performed at
		      	meta-feature level (eg. gene level).
	    
	    -O        	If specified, reads (or fragments if -p is specified) will
		      	be allowed to be assigned to more than one matched meta-
		      	feature (or feature if -f is specified). 
	    
	    -s <int>  	Indicate if strand-specific read counting should be performed.
		      	It has three possible values:  0 (unstranded), 1 (stranded) and
		      	2 (reversely stranded). 0 by default.
	    
	    -M        	If specified, multi-mapping reads/fragments will be counted (ie.
		      	a multi-mapping read will be counted up to N times if it has N
		      	reported mapping locations). The program uses the NH' tag to
		      	find multi-mapping reads.
	    
	    -Q <int>  	The minimum mapping quality score a read must satisfy in order
		      	to be counted. For paired-end reads, at least one end should
		      	satisfy this criteria. 0 by default.
	    
	    -T <int>  	Number of the threads. 1 by default.
	    
	    -R        	Output read counting result for each read/fragment. For each
		      	input read file, read counting results for reads/fragments will
		      	be saved to a tab-delimited file that contains four columns
		      	including read name, status(assigned or the reason if not
		      	assigned), name of target feature/meta-feature and number of
		      	hits if the read/fragment is counted multiple times. Name of
		      	the file is the same as name of the input read file except a
		      	suffix '.featureCounts' is added.
	    
	    --minReadOverlap <int>      Specify the minimum number of overlapped bases
		      	required to assign a read to a feature. 1 by default. Negative
		      	values are permitted, indicating a gap being allowed between a
		      	read and a feature.
	    
	    --readExtension5 <int>      Reads are extended upstream by <int> bases from
		      	their 5' end.
	    
	    --readExtension3 <int>      Reads are extended upstream by <int> bases from
		      	their 3' end.
	    
	    --read2pos <5:3>            The read is reduced to its 5' most base or 3'
		      	most base. Read summarization is then performed based on the
		      	single base which the read is reduced to.
	    
	    --fraction	If specified, a fractional count 1/n will be generated for each
		      	multi-mapping read, where n is the number of alignments (indica-
		      	ted by 'NH' tag) reported for the read. This option must be used
		      	together with the '-M' option.
	    
	    --primary 	If specified, only primary alignments will be counted. Primary
		      	and secondary alignments are identified using bit 0x100 in the
		      	Flag field of SAM/BAM files. All primary alignments in a dataset
		      	will be counted no matter they are from multi-mapping reads or
		      	not ('-M' is ignored). 
	    
	    --ignoreDup                 If specified, reads that were marked as
		      	duplicates will be ignored. Bit Ox400 in FLAG field of SAM/BAM
		      	file is used for identifying duplicate reads. In paired end
		      	data, the entire read pair will be ignored if at least one end
		      	is found to be a duplicate read.
	    
	    --countSplitAlignmentsOnly  If specified, only split alignments (CIGAR
		      	strings containing letter 'N') will be counted. All the other
		      	alignments will be ignored. An example of split alignments is
		      	the exon-spanning reads in RNA-seq data.
	    
	    Optional paired-end parameters:
	    
	    -p        	If specified, fragments (or templates) will be counted instead
		      	of reads. This option is only applicable for paired-end reads.
		      	The two reads from the same fragment must be adjacent to each
		      	other in the provided SAM/BAM file.
	    
	    -P        	If specified, paired-end distance will be checked when assigning
		      	fragments to meta-features or features. This option is only
		      	applicable when -p is specified. The distance thresholds should
		      	be specified using -d and -D options.
	    
	    -d <int>  	Minimum fragment/template length, 50 by default.
	    
	    -D <int>  	Maximum fragment/template length, 600 by default.
	    
	    -B        	If specified, only fragments that have both ends 
		      	successfully aligned will be considered for summarization.
		      	This option is only applicable for paired-end reads.
	    
	    -C        	If specified, the chimeric fragments (those fragments that 
		      	have their two ends aligned to different chromosomes) will
		      	NOT be included for summarization. This option is only 
		      	applicable for paired-end read data.
	    
	    -v        	Output version of the program.
	    
	    --donotsort   If specified, paired end reads will not be reordered even if
		      	reads from the same pair were found not to be next to each other
		      	in the input. 
	"> /dev/null
	echo "featurecounts done"
	########################################################################################################################
	date
}

function do_featurecounts_merge () {
	date
	########################################################################################################################
        sids="/data/NHLBI_BCB/Levine_Stew_Lab/06_Will_RNA-seq/sids.txt"
	out_dir="featurecounts"

	for i in `cat $sids | cut -f1 | head -1`; do
           SAMPLE=$i
	   echo -e "Gene\t$i" > tmp1
	   cat $SAMPLE/$out_dir/$i.genefeatureCounts.txt | sed '1,2d' | sort -k1,1 | cut -f 1,7 >> tmp1
	done

	for i in `cat $sids | cut -f1 | sed '1,1d' `; do
           SAMPLE=$i
	   echo -e "$i" > tmp2
	   cat $SAMPLE/$out_dir/$i.genefeatureCounts.txt | sed '1,2d' | sort -k1,1 | cut -f 7 >> tmp2
	   paste tmp1 tmp2 > tmp3
	   mv -f tmp3 tmp1
	done

	mv -f tmp1  gene.featurecount.txt
	rm -f tmp2
	rm -f tmp3
	########################################################################################################################
	date
}

function do_get_hisat2_stats () {
	date
	########################################################################################################################
	# change log files directory to absolute path of logfiles 
	logfiles="/data/NHLBI_BCB/Levine_Stew_Lab/06_Will_RNA-seq/02-fastqc-trimmomatic-hisat2-featurecounts/02-logfiles-trimmomatic-fastqc-hisat2-featurecounts"	

	cd $logfiles

	ls *.slurm.err.txt  | cut -f 1 -d'.' > tmp1
	grep 'reads; of these' *.slurm.err.txt | awk '{print $1}' | cut -f 2 -d':' | awk '{ printf("%'"'"'d\n",$1); }' > tmp2
	grep 'overall alignment rate' *.slurm.err.txt | awk '{print $1}' | cut -f 2 -d':' > tmp3
	grep 'aligned concordantly exactly 1 time' *.slurm.err.txt | awk '{print $2}' | awk '{ printf("%'"'"'d\n",$1); }' > tmp4
	grep 'aligned concordantly exactly 1 time' *.slurm.err.txt | awk '{print $3}' | sed 's/(//g' | sed 's/)//g' > tmp5
	echo -e "Sample\tTotal_Reads\tOverall_Alignment_Rate\tUniq_alignment\tUniq_alignment_%" > alignmetn.stats
	paste tmp1 tmp2 tmp3 tmp4 tmp5 >> alignmetn.stats
	rm -f tmp*
	a1=`cat alignmetn.stats | sed 's/,//g' | sed '1,1d'  | awk '{ sum += $2; n++ } END { print int(sum / n); }' | awk '{ printf("%'"'"'d\n",$1); }'`
	a2=`cat alignmetn.stats | sed 's/,//g' | sed '1,1d'  | awk '{ sum += $3; n++ } END { print sum / n; }' | awk -F. '{print $1"."substr($2,1,2)"%"}'` 
	a3=`cat alignmetn.stats | sed 's/,//g' | sed '1,1d'  | awk '{ sum += $4; n++ } END { print int(sum / n); }' | awk '{ printf("%'"'"'d\n",$1); }'`
	a4=`cat alignmetn.stats | sed 's/,//g' | sed '1,1d'  | awk '{ sum += $5; n++ } END { print sum / n; }' | awk -F. '{print $1"."substr($2,1,2)"%"}' `
	echo -e "Average\t$a1\t$a2\t$a3\t$a4" >> alignment.stats
	########################################################################################################################
	date
}

function do_get_featurecount_stats () {
	date
	########################################################################################################################
	sids="/data/NHLBI_BCB/Levine_Stew_Lab/06_Will_RNA-seq/sids.txt"
	
	# change log files directory to absolute path of logfiles 
	logfiles="/data/NHLBI_BCB/Levine_Stew_Lab/06_Will_RNA-seq/02-fastqc-trimmomatic-hisat2-featurecounts/02-logfiles-trimmomatic-fastqc-hisat2-featurecounts"	

	cd $logfiles
	
	echo -e "Sample\tGene_Count\tTotal_fragments\tSuccessfully_assigned\t%" > gene.featurecount.summary.txt

	for i in `cat $sids | cut -f1`; do
	    gene=`cat $i.slurm.err.txt | grep 'Meta-features : ' | awk '{print $4}'`
	   total=`cat $i.slurm.err.txt | grep 'Total fragments : ' | awk '{print $5}'`
	       s=`cat $i.slurm.err.txt | grep 'Successfully assigned fragments : ' | awk '{print $6"\t"$7}'`
	  echo -e "$i\t$gene\t$total\t$s" >> gene.featurecount.summary.txt
	done
	########################################################################################################################
	date
}

# do_fastqc
# do_trimmomatic
# do_fastqc
# do_hisat2
# do_featurecounts
# do_featurecounts_merge
# do_stringtie
do_get_hisat2_stats
# do_get_featurecount_stats




