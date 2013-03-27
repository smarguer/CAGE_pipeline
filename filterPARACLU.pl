#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
my $line1;
my $count=0;
my $count1=0;
my @holder;
my $chr=0;
my $start=0;
my $end=0;
my $strand="+";
my $first;
my $second;
my %read;
my %out;

my $length=500;
my $ratio=2;

#######################
#Takes input
#######################

if (@ARGV != 1) {die 'I need a paraclu output.';}
(my $in)=@ARGV;

open (IN, $in) or die 'could not find the paraclu output';
open (OUT, ">", $in.'.filtered') or die 'could not open output file.';

my $line=<IN>;
######################################
#parse paraclu output
######################################

while ($line=<IN>)
{
 $count++;
 chomp ($line);
 @holder = split (/\t/, $line);
################################################################################
 if(($holder[0] == $chr) && ($holder[1] eq $strand) && $holder[2] < $start) {die 'start values are not sorted correctly...';}
 
 if ((($holder[3]-$holder[2]) <= $length) && (($holder[7]/$holder[6]) >= $ratio))
 {
  if (($holder[0] != $chr) || ($holder[1] ne $strand) || ($holder[3] > $end)) 
  {
   print OUT "$line\n";
   $chr=$holder[0];
   $strand=$holder[1];
   $start=$holder[2]; 

   if($holder[1] eq $strand)
   {
    $end=$holder[3];
   }
   else
   {
    $end=0;  
   }
  }
 }
}

close IN;
close OUT;

