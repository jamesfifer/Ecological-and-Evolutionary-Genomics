# Ecological-and-Evolutionary-Genomics
# The following is the mapping and annotation pipeline for Symbiodiniaceae expression in "GOING WITH THE FLOW: CORALS IN HIGH-FLOW ENVIRONMENTS CAN BEAT THE HEAT"

# Reference assembly
1) Downloaded Cladocopium transcriptome from Ladner et al., 2012. https://www.ncbi.nlm.nih.gov/nuccore/GAFO01000001.1/

2) Ran ref.sh to recreate reference with Trinity headers via Trinity's genome guided assembly

3) Use the parallelBlast.sh to blast against NCBI nt database and retains hits that match Dinophyceae with an evalue cutoff of 1e-5. 
Filter out contigs with <300bp using the 300bpfilter.sh & contigs_to_fasta.py

4) Use mapping.sh to generate counts file via RSEM. Take "raw counts" from RSEM gene file. 

# Annotation 
1) Use the pipeline outlined in Annotation.sh

# DESeq and analysis
1) Use SymFlow.rmd file to run deseq2 and downstream analyses

# GO MWU
1) Use the GO_MWU.rmd file taken from https://github.com/z0on/GO_MWU to generate GO enrichment plots




