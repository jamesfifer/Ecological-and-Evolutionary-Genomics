#first filter the blastoutput file so it only contains bp >300
awk '$4 >300 '{print $1} > blastoutput.300bp.plus.txt
#then use grep -F -f to filter out the IDs that are shared between
grep -F -f blastoutput.300bp.plus.txt SymContigs.txt > filteredSymContigs.txt
#then put the filteredtxt file through the Contigs_to_fasta.py script 

