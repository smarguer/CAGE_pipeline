#!/usr/bin/perl

use strict;
use warnings;

my $line;
my $count=0;
my @holder;
my $gene;
my %out;

if (@ARGV != 1) {die "wrong number of files";}
(my $in)=@ARGV;

open (IN, $in) or die 'could not find the input file';

while ($line=<IN>)
 {
  $count++;
  chomp ($line);
  
  @holder = split (/\t/, $line);
  
  $gene="$holder[0]\t$holder[1]\t$holder[2]\t$holder[3]\t$holder[4]\t$holder[5]";
  if(defined($out{$gene}{CL}))
  {
   $out{$gene}{CL}=$holder[6];
   $out{$gene}{CH}=$holder[7];
   $out{$gene}{ST}=$holder[8];
   $out{$gene}{S}=$holder[9];
   $out{$gene}{E}=$holder[10];
   $out{$gene}{CL2}=$holder[11];
   $out{$gene}{N}=$holder[12];
   $out{$gene}{M}=$holder[13];
  }
  else
  {
   $out{$gene}{CL}.=";$holder[6]";
   $out{$gene}{CH}.=";$holder[7]";
   $out{$gene}{ST}.=";$holder[8]";
   $out{$gene}{S}.=";$holder[9]";
   $out{$gene}{E}.=";$holder[10]";
   $out{$gene}{CL2}.=";$holder[11]";
   $out{$gene}{N}.=";$holder[12]";
   $out{$gene}{M}.=";$holder[13]";
  }
 
 }
 for my $out (keys %out)
 {
  print "$out\t$out{$gene}{CL}\t$out{$gene}{CH}\t$out{$gene}{ST}\t$out{$gene}{S}\t$out{$gene}{E}\t$out{$gene}{CL2}\t$out{$gene}{N}\t$out{$gene}{M}\n";
 }

close IN;

