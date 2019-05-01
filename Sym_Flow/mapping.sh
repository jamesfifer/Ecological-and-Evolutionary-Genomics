#!/bin/bash

#below code most be run while in directory with samples  otherwise you will get unable to locate sample__
#Trinity recommends running align and estimate abundance first without sample file to prep reference and then with sample file without -prep reference command. (Prep the reference and then run on samples only). 

module load bowtie2/2.3.2
module load rsem/1.3.0
module load samtools/1.9


#/projectnb/bi594/jfifer/final_project/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl --transcripts /projectnb/bi594/jfifer/final_project/trinity_out_dir/Trinity-GG.fasta --est_method RSEM --aln_method bowtie2 --prep_reference --trinity_mode

/projectnb/bi594/jfifer/final_project/trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl --transcripts /projectnb/bi594/jfifer/final_project/trinity_out_dir/Trinity-GG.fasta --est_method RSEM --aln_method bowtie2  --trinity_mode --output_dir ./DE_output --samples_file /projectnb/bi594/jfifer/final_project/samples.txt --seqType fq




