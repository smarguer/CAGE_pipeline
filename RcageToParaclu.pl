#!/usr/bin/perl

use strict;
use warnings;


if (@ARGV != 1) {die "wrong number of files";}
(my $in)=@ARGV;

open (IN, $in) or die 'could not find the input file';

my $line=<IN>;
my $count=0;
my @holder;
my %chr=("chrI" => 1, "chrII" => 2,"chrIII" => 3, "chrMT" => 4, "chrMTR" => 5, "chrAB325691" => 6);

while ($line=<IN>)
 {
  $count++;
  chomp ($line);

  @holder = split (/\t/, $line);
  print "$chr{$holder[1]}\t$holder[4]\t$holder[2]\t$holder[3]\t$holder[0]\t$holder[5]\t$holder[6]\n"  

 }


close IN;

