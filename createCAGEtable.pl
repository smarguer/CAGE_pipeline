#!/usr/bin/perl

use strict;
use warnings;
use POSIX;
my $line;
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
my %clusters;
my %out;
my @first=();
my @second=();
my $match=0;
my $antisense;
my $name;
my $add_name;
#######################
#Takes input
#######################

if (@ARGV != 2) {die 'I need a paraclu mapped file  and a gff file';}
(my $in,my $gff)=@ARGV;

open (IN, $in) or die 'could not find the paraclu output';
open (GFF, $gff) or die 'could not find the gff file';
my $line1=<GFF>;

######################################
#parses Exonerate output
######################################

while ($line=<IN>)
{
 chomp ($line);
 @holder = split (/\t/, $line);
 
 unless(exists $clusters{$holder[10]})
 {
  $clusters{$holder[10]}=$line;
 }
 else
 {
  $count++;
  $add_name=$holder[10].$count;
  $clusters{$add_name}=$line;
 }
}

while ($line1=<GFF>)
{
 $count++;
 chomp ($line1);
 @holder1 = split (/\t/, $line1);
 if ($holder1[2] !~ /gene/)
 {
  next;
 }
 if ($holder1[9] =~ /LTR/)
 {
  $name = 'LTR'.'.'.$holder1[3].'.'.$holder1[12];
 }
 else
 {
  $name=$holder1[9];
 }
 $antisense='AS.'.$name;

################################################################################
 foreach my $out (keys %clusters)
 {
  #print "$out\t###\t$holder1[9]\n";
  if (($out =~ /$name\./) && ($out !~ /$antisense/))
  {
   print "$name\t$holder1[3]\t$holder1[4]\t$holder1[6]\t$holder1[12]\t$holder1[8]\t$clusters{$out}\n";
   $match=1;
  }
 }
 if ($match==0)
 {
   print "$name\t$holder1[3]\t$holder1[4]\t$holder1[6]\t$holder1[12]\t$holder1[8]\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
 }
 $match=0;
}

close IN;
close GFF;
