# MSAvariants
This script takes a multiple sequence alignment in fasta format as input and identifies all single nucleotide variants present in the alignment. Missing and/or poor quality bases can be encoded as "N" and will be ignored. This script requires that BioPerl be installed. 
USAGE: perl identify_MSA_variants.pl ALIGNMENT.FASTA > VARIANTS.OUT
