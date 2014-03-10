#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
my $line;
my $line1;
my $count=0;
my $count1=0;
my @holder;
my @holder1;
my $chr;
my $start;
my $end;
my $strand;
my $first;
my $second;
my %read;
my %out;
my @first=();
my @second=();

#######################
#Takes input
#######################

if (@ARGV != 2) {die 'I need a paraclu output and a reads "fused" file';}
(my $in,my $in1)=@ARGV;

open (IN, $in) or die 'could not find the paraclu output';
open (IN1, $in1) or die 'could not find the fused file output';

######################################
#parses Exonerate output
######################################

while ($line1=<IN1>)
{
 chomp ($line1);
 @holder1 = split (/\t/, $line1);
 $read{$holder1[8]}{F}=$holder1[7];
 $read{$holder1[8]}{S}=$holder1[6];
 $read{$holder1[8]}{CH}=$holder1[10];
 $read{$holder1[8]}{ST}=$holder1[9];
}

while ($line=<IN>)
{
 $count++;
 chomp ($line);
 @holder = split (/\t/, $line);
################################################################################
 $chr=$holder[0];
 $strand=$holder[1];
 $start=$holder[2]; 
 $end=$holder[3]; 

 foreach my $out (keys %read)
 {
  if (($read{$out}{CH}==$chr)&&($read{$out}{ST} eq $strand)&&($out >= $start)&&($out <= $end))
  {
   push(@first,$read{$out}{F});
   push(@second,$read{$out}{S});
  }   
 }
 %out = map { $_, 1 } @first;
 $first=join(';',keys %out);
 %out=();
 %out = map { $_, 1 } split (';',$first);
 $first=join(';',keys %out);
 %out=();

 %out = map { $_, 1 } @second;
 $second=join(';',keys %out);
 %out=();
 %out = map { $_, 1 } split (';',$second);
 $second=join(';',keys %out);
 %out=();

 unless(defined $first && length $first)
 {
  $first="OUTSIDE";
 }
 unless(defined $second && length $second)
 {
  $second="OUTSIDE";
 }
 

 %out=();
 @first=();
 @second=();

 print "c".$count."\t$line\t$first\t$second\n";
}

close IN;
close IN1;
