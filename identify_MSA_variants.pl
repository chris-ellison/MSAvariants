#This script takes a multiple sequence alignment in fasta format as input and identifies all single nucleotide variants present in the alignment. Missing and/or poor quality bases can be encoded as "N" and will be ignored. This script requires that BioPerl be installed. 
#USAGE: perl identify_MSA_variants.pl ALIGNMENT.FASTA > VARIANTS.OUT

use strict;
use warnings;
use Bio::SeqIO;

my $msa = shift;
my @seq_matrix;
my %seq_length;
my $seq_count=0;

print "\#POSITION\tALLELES";
my $seqin = Bio::SeqIO->new( -format => 'FASTA' , -file => $msa);
while((my $seqobj = $seqin->next_seq())) {
    my $string = $seqobj->seq();
    my $seq_id = $seqobj->display_id();
    $string =~ tr/atcgn/ATCGN/;
    $seq_id =~ s/\s+/_/g;
    my @columns = split("",$string);
    my $arLength = scalar(@columns);
    $seq_length{$arLength}++;
    push(@seq_matrix,[@columns]);
    print "\t",$seq_id;
    $seq_count++;
}
print "\n";

my @each_length = keys %seq_length;
die "SEQUENCES IN ALIGNMENT ARE DIFFERENT LENGTHS!\n" if scalar(@each_length)>1;
my $len = $each_length[0];

for(my $col=0;$col<$len;$col++){
    my %alleles;
    my $col_string;
    for(my $row=0;$row<$seq_count;$row++){
	my $genotype = $seq_matrix[$row][$col];
	if($genotype ne "N"){
	    $alleles{$genotype}++;
	}
	$col_string.="\t".$genotype;
    }
    my @all = sort {$alleles{$b}<=>$alleles{$a}} keys %alleles;
    next unless scalar(@all)>1;
    print $col+1,"\t",$all[0],":",$alleles{$all[0]};
    for(my $i=1;$i<scalar(@all);$i++){
	print ",",$all[$i],":",$alleles{$all[$i]};
    }
    print $col_string,"\n";
}
