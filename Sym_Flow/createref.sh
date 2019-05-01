#!/bin/bash -l
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N createref # job name, anything you want
#$ -pe omp 4
#$ -l h_rt=44:00:00 #maximum run time
#$ -P bi594
#$ -M james.e.fifer@gmail.com #your email
#$ -m be

module load boost/1.58.0
module load bowtie2
module load tophat/2.1.1
module load samtools/1.9

#use bowtie2 to build the database before running this script



#for sample in `ls ./*R1_Combined.fq`
#do
#	outdir=`echo ${sample/_R1_Combined.fq/}`
#	tophat -p 4 -o $outdir /projectnb/bi594/jfifer/final_project/Cladocopium_genome $sample ${sample/R1/R2}
#done

#add -p for more threads
#Takes all files in all downhill directories with name and copies with unique names to current directory

#find ./ -name 'accepted_hits.bam' -exec cp --backup=numbered -t ./ {} + 
#samtools merge merged.bam *bam*

#Trinity recommends to make sure the file is coordinate sorted by running 'samtools sort' on it

#samtools sort merged.bam > accepted_hits.merged.sorted.bam

#Run genome guided trinity assembly
module load java/1.8.0_151
module load trinity/2.4.0
Trinity --genome_guided_bam /projectnb/bi594/jfifer/final_project/accepted_hits.merged.sorted.bam \
       --CPU 2 --max_memory 1G --genome_guided_max_intron 5000
