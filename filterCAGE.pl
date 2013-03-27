#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
my $line;
my $count=0;
my $count1=0;
my @holder;
my $limit=1000;
my $read;
my %read;
my %read1;
#######################
#Takes input
#######################

if (@ARGV != 1) {die 'I need reads "fused" file';}
(my $in)=@ARGV;

open (IN, $in) or die 'could not find the paraclu output';
open (OUT, ">", $in.'.filtered') or die 'could not open output file';
open (OUT1, ">", $in.'.filtered.paracluIN') or die 'could not open output file';
######################################
#parses Exonerate output
######################################

while ($line=<IN>)
{
 chomp ($line);
 @holder = split (/\t/, $line);
##removes pairs with single read mapping
 if ($holder[12] =~ /sin/)
 {
  next;
 }
##removes pairs within multiple matches
 if (($holder[5] eq "M")||($holder[11] eq "M"))
 {
  next;
 }
##removes rRNA
 if ($holder[6] =~ /SPRRNA/)
 {
  next;
 }
 if ($holder[7] =~ /SPRRNA/)
 {
  next;
 }
##removes chimeras, adjust $limit
 if ($holder[13] > $limit)
 {
  next;
 }
##removes replicated pairs
 if($holder[3] eq "+")
 {
  $read=$holder[1].$holder[8];
 }
 elsif ($holder[3] eq '-')
 {
  $read=$holder[2].$holder[8];
 }
 if(exists $read{$read})
 {
  next;
 }
 else
 {
  $read{$read}=1;
  print OUT "$line\n";
 }
##creates hash for paraclu input file

 if(exists $read1{$holder[10]}{$holder[9]}{$holder[8]})
 {
  $read1{$holder[10]}{$holder[9]}{$holder[8]}++;
 }
 else
 {
  $read1{$holder[10]}{$holder[9]}{$holder[8]}=1;
 }
}

##prints out paraclu input file
for my $c_out (sort { lc($a) cmp lc($b) } keys %read1)
{
 for my $s_out (keys %{$read1{$c_out}})
 {
  for my $start_out (sort { $a <=> $b } keys %{$read1{$c_out}{$s_out}})
  {
   print OUT1 "$c_out\t$s_out\t$start_out\t$read1{$c_out}{$s_out}{$start_out}\n";
  }
 }
}


close IN;
close OUT;
close OUT1;















