#!/usr/bin/perl


use strict;
use warnings;
use POSIX;

#######################
#Takes input
#######################

if (@ARGV != 3) {die "I need first the read 1 and then the read 2 mapping outputs";}
#(my $in1,my $in2)=@ARGV;
(my $in1,my $in2,my $in3)=@ARGV;

open (IN1, $in1) or die 'could not find the read 1 mapping output';
open (IN2, $in2) or die 'could not find the read 2 read annotation';
open (IN3, $in3) or die 'could not find the read 1 read annotation';

$in1 =~ s{.*/}{};
open (OUT, ">", './'.$in1.'.fused') or die 'could not open output file...';

#######################
#variables
#######################
my $line1;
my $line2=<IN2>;
my $line3=<IN3>;
my $count=0;
my $count1=0;
my @holder;
my $names;
my $chr;
my $start;
my $strand;
my $read;
my $unique;
my $flag;
my $insert;
my $mid;
my $map;
my %read1;
my %read2;
my @read;
my @output1;
my @output2;
my @output3;

########################
#reads read 2 annotation
########################
while ($line2=<IN2>)
{
	
 chomp ($line2);
 @holder = split (/\t/, $line2);
 unless ($holder[0]){next;}
 $names=join(';',@holder[5..$#holder]);
#print "$names\n";
 @read=split /;/, $holder[0];    
 $read=substr $read[0],1;
 $unique=substr $holder[0],0,1;
#print "$read\n";
 $read2{$read}="$read\t".join("\t",@holder[1..4])."\t$unique\t$names";
#print "$read{$read}\n";

}
########################
#reads read 1 annotation
########################
while ($line3=<IN3>)
{
	
 chomp ($line3);
 @holder = split (/\t/, $line3);
 unless ($holder[0]){next;}
 $names=join(';',@holder[5..$#holder]);
#print "$names\n";
 @read=split /;/, $holder[0];    
 $read=substr $read[0],1;
 $unique=substr $holder[0],0,1;
#print "$read\n";
 $read1{$read}="$read\t".join("\t",@holder[1..4])."\t$unique\t$names";
#print "$read{$read}\n";

}


#######################################
#parses Exonerate output
######################################

###############################################################################
while ($line1=<IN1>)
{
 chomp ($line1);     
 @holder = split (/ /, $line1);
 unless ($holder[0]){next;}
 if ($holder[0] !~ /vulgar/)
 {next;}
 @read=split /;/, $holder[1];    
 $holder[5] =~ /CH(\d{1})/;
 $chr=$1;
 $unique=substr $holder[0],0,1;
 if(exists $read2{$read[0]})
 {
  @output2=split (/\t/, $read2{$read[0]});
  if($output2[3] eq "+")
  {
   $mid=$output2[1];
  }
  elsif($output2[3] eq "-")
  {
   $mid=$output2[2];
  }
  if($output2[0] eq $read[0])
  {
   $flag="ok";
  }
  else
  {
   $flag="dif";
  }
 }
 else
 {
  @output2=($read[0],"NA","NA","NA","NA","NA","NA");
  $mid=0;
  $flag="sin";
 }
 $strand=$holder[8];

 if ($strand eq "+")
 {
  $start=$holder[6]; 
 }
 elsif ($strand eq "-")
 {
  $start=$holder[7]; 
 }    
 if($mid != 0)
 {
  $insert=abs($start-$mid);
 }
 else
 {
  $insert="NA";
 }
 @output1=($start,$strand,$chr,$unique,$flag,$insert);
 @output3=split (/\t/, $read1{$read[0]});
 if($output2[6] eq $output3[6])
 {$map="same";}
 else
 {$map="dif";}
 

 push(@output2,$output3[6]);
 push(@output2,@output1);
 push(@output2,$map);
 print OUT join("\t",@output2)."\n";
}
close IN1;
close IN2;
close IN3;
close OUT;
