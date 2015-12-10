#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Bio::PopGen::Statistics;
use Bio::AlignIO;
use Bio::PopGen::IO;
use Bio::PopGen::Simulation::Coalescent;

my($input,$sites,$out,$list);
GetOptions(
  "input|i:s"  =>  \$input,
  "number|n:i" =>  \$sites,
  "out|o:s"    =>  \$out,
  "list|l:s"   =>  \$list,
  "help|h"     =>  sub{&usage;exit(-1);}
);

unless(defined $input && defined $sites){&usage;exit(-1);}
my $OUT;
unless(defined $out){$OUT = \*STDOUT;}
else{open $OUT,">$out";}

# Only for individuals in --list
my @pop;
if (defined $list){
  open FH,$list;
  while(<FH>){
    chomp;
        push @pop,$_;
  }
}

my $stats = Bio::PopGen::Statistics->new();
my $io = Bio::PopGen::IO->new(-format => 'prettybase',-file  => $input);

print $OUT "#Name\tPi\tTheta\tTajima_D\tFu_and_Li_D*\tFu_and_Li_F*\n";
print $OUT "$input\t";

unless (@pop>=1){
  if(my $pop = $io->next_population) {
    my $PI = $stats->pi($pop,$sites);
    my $THETA  = $stats->theta($pop,$sites);
    my $tajima_D = $stats->tajima_D($pop);
    my $fu_and_li_D_star = $stats->fu_and_li_D_star($pop);
    my $fu_and_li_F_star = $stats->fu_and_li_F_star($pop);
        print $OUT "$PI\t$THETA\t$tajima_D\t$fu_and_li_D_star\t$fu_and_li_F_star\n";
  }
}else{
  my $PI = $stats->pi(\@pop,$sites);
  my $THETA  = $stats->theta(\@pop,$sites);
  my $tajima_D = $stats->tajima_D(\@pop);
  my $fu_and_li_D_star = $stats->fu_and_li_D_star(\@pop);
  my $fu_and_li_F_star = $stats->fu_and_li_F_star(\@pop);
  print $OUT "$PI\t$THETA\t$tajima_D\t$fu_and_li_D_star\t$fu_and_li_F_star\n";
}

sub usage{
print STDERR <<HELP
Usage:
  perl $0 -i [file.prettybase] -n [number of sites] -l <indis.list> -o <out>
Parats:
  --input|-i    file in prettybase format
  --number|-n   number of nucleotide sites or the sequence length
  --list|-l     file with each row corresponding an individul
  --out|-o      output file name
Funct:
  Return the pi(per site), theta_W(per site), Tajima's D, Fu and Li's D*/F*
  for all individuals or only individuals in --list
HELP
}
