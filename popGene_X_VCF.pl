#!/usr/bin/perl
use warnings;
use strict;
use List::Util qw(min);
use Getopt::Long;
use Bio::PopGen::Statistics;
use Bio::AlignIO;
use Bio::PopGen::IO;
use Bio::PopGen::Simulation::Coalescent;

my ($vcf,$chrsize,$gender,$window);
$window = 3000000;
GetOptions(
  "vcf|v:s"    =>    \$vcf,
  "chr|c:s"    =>    \$chrsize,
  "win|w:i"    =>    \$window,
  "gender|g:s" =>    \$gender,
  "help|h"     =>    sub{&usage;exit(-1);}
);

unless(defined $chrsize && defined $vcf && defined $gender){&usage;exit(-1);}
else{open CHR,$chrsize; open VCF,$vcf; open GEN,$gender;}

my %chr;
while(<CHR>){
  chomp;
  my @t = split;
  $chr{$t[0]} = $t[1];
}

my %gender;
while(<GEN>){
  chomp;
  my @t = split;
  $gender{"$t[0]"} = $t[1];
}

my @invi;
my ($chr,$start,$end);
my @cSW;

while(<VCF>){
  chomp;
  my $line = $_;
  next if /^##/;
  my @t = split;
  if(/^#/){
    @invi = @t[9..$#t];
    next;
  }

core:{

  unless(defined $chr){
    $chr = $t[0];
    $start = 0;
    $end = min($window,$chr{$chr});
    goto core;
  }else{
    if($chr eq $t[0]){
      if($t[1]<=$start){
        print STDERR "Errors: non-sorted vcf file.\n";
        exit(-1);
      }elsif($t[1]>$start && $t[1]<=$end){
        push @cSW,$line;
      }else{
        if (@cSW == 0){
          print "$chr\t$start\t$end\t0\t0\t0\t0\t0\n";
        }else{
          &vcf2prettybase(@cSW);
          &popGene(join("_",($chr,$start,$end)),$end-$start);
        }
        @cSW = ();
        ($chr,$start,$end) = &nextregion($chr,$start,$end);
        goto core;
      }
   }else{
     if (@cSW == 0){
       print "$chr\t$start\t$end\t0\t0\t0\t0\t0\n";
     }else{
       &vcf2prettybase(@cSW);
       &popGene(join("_",($chr,$start,$end)),$end-$start);
     }
     @cSW = ();
# report the regions (with no polymorphisms) near to the end of the chromosome if exists
     other:{
       ($chr,$start,$end) = &nextregion($chr,$start,$end);
       if(defined $chr){
         print "$chr\t$start\t$end\t0\t0\t0\t0\t0\n";
         goto other;
       }
     }
     goto core;
    }
  }
} # core

}


if (@cSW == 0){
  print "$chr\t$start\t$end\t0\t0\t0\t0\t0\n";
}else{
  &vcf2prettybase(@cSW);
  my $file =  join("_",($chr,$start,$end));
  &popGene($file,$end-$start);
}

# report the regions (with no polymorphisms) near to the end of the chromosome if exists
other2:{
  ($chr,$start,$end) = &nextregion($chr,$start,$end);
  if(defined $chr){
    print "$chr\t$start\t$end\t0\t0\t0\t0\t0\n";
    goto other2;
  }
}

## subrountines

sub nextregion{
  my ($C,$S,$E) = @_;
  if ($E != $chr{$C}){
    return ($chr,$E,min($E+$window,$chr{$C}));
  }else{
    return (undef,undef,undef);
  }
}

sub vcf2prettybase{

  my $out = join("_",($chr,$start,$end));
  open OUT, ">$out";

  foreach(my $i=0;$i<=$#_;$i++){
    my @t = split("\t",$_[$i]);
    foreach(my $j=9;$j<=$#t;$j++){
      my $g = $invi[$j-9];
          if($gender{$g} =~ /^M/i){
        if($t[$j]=~/1/){
          print OUT "A","$i","\t$g","\t$t[4]\n";
        }else{
          print OUT "A","$i","\t$g","\t$t[3]\n";
        }
      }elsif($gender{$g} =~ /^F/i){
        my @geno = ($1,$3) if $t[$j] =~ /^(.)(.)(.)/;
        $geno[0] = (($geno[0] eq 1) ? $t[4] : $t[3]);
        $geno[1] = (($geno[1] eq 1) ? $t[4] : $t[3]);
        print OUT "A","$i","\t$g",".1","\t$geno[0]\n";
        print OUT "A","$i","\t$g",".2","\t$geno[1]\n";
      }
    }
  }
  close OUT;
}


sub popGene{
  my ($input,$sites) = @_;
  my $stats = Bio::PopGen::Statistics->new();
  my $io = Bio::PopGen::IO->new(-format => 'prettybase',-file  => $input);

  if(my $pop = $io->next_population) {
    my $PI = $stats->pi($pop,$sites);
    my $THETA  = $stats->theta($pop,$sites);
    my $tajima_D = $stats->tajima_D($pop);
    my $fu_and_li_D_star = $stats->fu_and_li_D_star($pop);
    my $fu_and_li_F_star = $stats->fu_and_li_F_star($pop);

    print "$chr\t$start\t$end\t$PI\t$THETA\t$tajima_D\t$fu_and_li_D_star\t$fu_and_li_F_star\n";
    system("rm $input");
  }
}

sub usage{
print STDERR <<HELP
Usage:
  perl $0 -v [vcf] -c [chr.size] -w [3000000]
Parats:
  --vcf|-v         please give ordered/pahsed vcf file for X chromosome
  --chr|-c         chr.size
  --win|-w         window size [3000000]
  --gender|-g      file containing the gender of each individual in vcf (ID\t[Mm](ale)*/{Ff](emale)*)
  --help|-h        print help message
Funct:
  Return Pi(Per Site),Theta(Per Site),Tajima's D, Fu and Li's D*/F* for each window
HELP
}
