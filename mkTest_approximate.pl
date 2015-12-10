#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(min);
use Bio::PopGen::Statistics;
use Bio::Tools::CodonTable;

unless(@ARGV==1){
  &usage;exit(-1);
}else{
  open FH,$ARGV[0];
}

# standard codon table
my $codonTable = Bio::Tools::CodonTable->new();
$codonTable->id(1);

my @counts = (0,0,0,0); # NP,NF,SP,SF

my @seq; # array of array
while(<FH>){
  next if /^>/;
  chomp;
  push @seq,uc($_);
}

for(my $i=0; $i<= length($seq[0])-3; $i+=3){
  my @p1 = (substr($seq[0],$i,1),
            substr($seq[1],$i,1),
            substr($seq[2],$i,1),
            substr($seq[3],$i,1));
  my @p2 = (substr($seq[0],$i+1,1),
            substr($seq[1],$i+1,1),
            substr($seq[2],$i+1,1),
            substr($seq[3],$i+1,1));
  my @p3 = (substr($seq[0],$i+2,1),
            substr($seq[1],$i+2,1),
            substr($seq[2],$i+2,1),
            substr($seq[3],$i+2,1));
  my @codon = (substr($seq[0],$i,3),
               substr($seq[1],$i,3),
               substr($seq[2],$i,3),
               substr($seq[3],$i,3));

  if(&uniq(@p1)+&uniq(@p2)+&uniq(@p3) == 4){

        next if $codonTable->translate($codon[0]) eq "*";
        next if $codonTable->translate($codon[1]) eq "*";
        next if $codonTable->translate($codon[2]) eq "*";
        next if $codonTable->translate($codon[3]) eq "*";

        if($codon[0] eq $codon[1] && $codon[2] eq $codon[3]){
      if($codonTable->translate($codon[0]) eq $codonTable->translate($codon[2])){
        $counts[3]++;
      }else{
        $counts[1]++;
      }
    }elsif($codon[0] eq $codon[1]){
      if($codonTable->translate($codon[2]) eq $codonTable->translate($codon[3])){
        $counts[2]++;
      }else{
        $counts[0]++;
      }
    }elsif($codon[2] eq $codon[3]){
      if($codonTable->translate($codon[0]) eq $codonTable->translate($codon[1])){
        $counts[2]++;
      }else{
        $counts[0]++;
      }
    }
  }
}

print "$ARGV[0]\t",join("\t",@counts),"\t";

if($counts[1]!=0 && $counts[2]!=0 && $counts[3]!=0){
  my $NI = ($counts[0]/$counts[2])/($counts[1]/$counts[3]);
  my $stats = Bio::PopGen::Statistics->new();
  print "$NI\t";
  print $stats->mcdonald_kreitman_counts(@counts);
  print "\n";
}else{
  print "NA\tNA\n";
}

sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}


sub usage{
print STDERR <<HELP
Usage:
  perl $0 seq.aln.
Funct:
  Approximate Version of the McDonald-Kreitman Test (similar to R package PopGenome).
Warns:
  0). Four sequences in seq.aln, including ref-spe1,collapsed-alt-spe1,ref-spe2 and collapsed-alt-spe2
  1). Codons with indels are already excluded.
  2). Stop codons are not counted
  3). Only Codons with single SNP are included (similar to R package PopGenome).
HELP
}
