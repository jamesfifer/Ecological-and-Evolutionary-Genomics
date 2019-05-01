##This pipeline was developed by Dr. Sarah Davies at Boston University
#######################################
Annotating the transcriptome

# getting uniprot_swissprot KB database
echo "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" >getuni2
nano getuni2
#copy and paste this text into the top of the file
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N getuni # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M daviessw@gmail.com #your email
#$ -m be
qsub getuni2

# getting annotations (this file is over 3G, will take a while)
echo "wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz" >getgo
nano getgo
#copy and paste this text into the top of the file
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N getgo # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M daviessw@gmail.com #your email
#$ -m be
qsub getgo

# unzipping
gunzip uniprot_sprot.fasta.gz &
gunzip idmapping_selected.tab.gz &

# indexing the fasta database
module load blast+/2.2.29
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
nano mbd
#copy and paste this text into the top of the file
#!/bin/bash
#$ -V # inherit the submission environment
#$ -cwd # start job in submission directory
#$ -N mbd # job name, anything you want
#$ -l h_rt=24:00:00 #maximum run time
#$ -M daviessw@gmail.com #your email
#$ -m be
qsub mbd

# splitting the transcriptome into 100 chunks
splitFasta.pl crepidula.fasta 100

# blasting all 40 chunks to uniprot in parallel, 3 cores per chunk
module load blast+/2.2.29
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 3 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
module load python3
./scc6_qsub_launcher.py -N blast -P bi594 -M daviessw@gmail.com -j y -h_rt 24:00:00 -jobsfile bl
qsub blast_array.qsub

# combining all blast results
cat subset*br > crepidula.br
rm -f subset*

# for trinity-assembled transcriptomes: annotating with "isogroup" (=component)
grep ">" crepidula.fasta | perl -pe 's/>((\S+)_i\d+).+/$1\t$2/' >crepidula_seq2iso.tab 
cat crepidula.fasta | perl -pe 's/>((\S+)_i\d+)/>$1 gene=$2/' >crepidula_iso.fasta

# extracting gene names (per isogroup):
echo "getGeneNameFromUniProtKB.pl blast=crepidula.br prefix=crepidula fastaQuery=crepidula_iso.fasta" >getgn
module load python3
./scc6_qsub_launcher.py -N getgn -P bi594 -M daviessw@gmail.com -j y -h_rt 24:00:00 -jobsfile getgn
qsub getgn_array.qsub

# extracting GO annotations (per isogroup)
echo "getGOfromUniProtKB.pl blast=crepidula.br prefix=crepidula fastaQuery=crepidula_iso.fasta" >getgo
module load python3
./scc6_qsub_launcher.py -N getgo -P bi594 -M daviessw@gmail.com -j y -h_rt 24:00:00 -jobsfile getgo
qsub getgo_array.qsub
