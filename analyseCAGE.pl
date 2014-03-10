use strict;
use warnings;
use Benchmark;
use Data::Dumper;
use POSIX qw(strftime);

my $scriptp="/jurg/homedir/samuel/POMBE_SEQ/analysis/SCRIPTS/CAGE_pipeline/";
my %status = ();
my $n = 1;
my @command=@ARGV;

my $use = "\nUse:\n-gp: gff path (/jurg/homedir/samuel/POMBE_SEQ/analysis/)\n-qp: mapped reads files path (default current directory)\n-g1: gff used for mapping (default gff_090511.txt)\n-g2: gff to be used for UTR call (default gff_090511.txt)\n-r1: read 1 file name (no default)\n-r2: read 2 file name (no default)\n-m: paraclu lower cluster signal limit (default = 30)\n\n";

##-gp
my $gffp = "/jurg/group/SAM_RNA-SEQ_PIPELINE/ANNOT/";
##-qp
my $fastqp = "./";
##-g
my $gff1 = "gff_011112.txt";
my $gff2 = "gff_090511.txt";
##-r1
my $read1;
##-r2
my $read2;
##-m
my $minClust=30;

##-t
#my $tag = $today[2].$transTime{$today[1]}.substr($today[4],2,2);

my $readinput;

die $use unless (defined(@ARGV));

while (defined($readinput=shift))
{

 die $use unless ($readinput =~ /-/);

 if($readinput eq '-gp')
 {
  $gffp=shift;
 }

 if($readinput eq '-qp')
 {
  $fastqp=shift;
 }

 if($readinput eq '-r1')
 {
  $read1=shift;
 }

 if($readinput eq '-r2')
 {
  $read2=shift;
 }

 if($readinput eq '-g1')
 {
  $gff1=shift;
 }

 if($readinput eq '-g2')
 {
  $gff2=shift;
 }

 if($readinput eq '-m')
 {
  $minClust=shift;
 }

# if($readinput eq '-t')
# {
#  $tag=shift;
# }
}

die "\ngff1 file with unexpected name structure please use ".'"'."gff_DDMMYY.txt".'"'."\n$use\n" unless ($gff1 =~ /gff_\d{6}\.txt/);
die "\ngff2 file with unexpected name structure please use ".'"'."gff_DDMMYY.txt".'"'."\n$use\n" unless ($gff2 =~ /gff_\d{6}\.txt/);
die "\n-r1 option missing with no default\n$use\n"  unless (defined($read1));
die "\n-r2 option missing with no default\n$use\n"  unless (defined($read2));
#print "$options\n$gffp\n$fastqp\n$fastq\n$gff\n$read\n$tag\n";
#print "Let's go!\n";

#####################################################################
#time it
my $time0 = new Benchmark;

my $now = localtime;

print "\nCommand: $0 @command\n\ndate: $now\n\n";


my $debug=0;

if ($debug==0)
{
print "fusing read 1 and read 2 files\n";
$status{$n}=system "perl ".$scriptp."Map1rtWith2nd.pl ".$fastqp.$read1.' '.$fastqp.$read2.'.mapIV.'.$gff1.' '.$fastqp.$read1.'.mapIV.'.$gff1;
die "Problem stage $n" if $status{$n}; $n++;

print "filtering...\n";
$status{$n}=system "perl ".$scriptp."filterCAGE.pl ".$read1.'.fused';
die "Problem stage $n" if $status{$n}; $n++;

print "running paraclu\n";
$status{$n}=system '/jurg/homedir/samuel/POMBE_SEQ/analysis/SCRIPTS/CAGE_pipeline/paraclu-9/paraclu '.$minClust.' '.$read1.'.fused.filtered.paracluIN > '.$read1.'.fused.filtered.paracluOUT';
die "Problem stage $n" if $status{$n}; $n++;

print "filtering paraclu output (removing clusters of length > 500nt, with density ratio < 2x, or included in larger ones\n";
$status{$n}=system "perl ".$scriptp."filterPARACLU.pl ".$read1.'.fused.filtered.paracluOUT';
die "Problem stage $n" if $status{$n}; $n++;

print "mapping paraclu clusters to annotation\n";
$status{$n}=system "perl ".$scriptp."mapCAGEclusters.pl ".$read1.'.fused.filtered.paracluOUT.filtered '.$read1.'.fused.filtered > '.$read1.'.CAGEmap';
die "Problem stage $n" if $status{$n}; $n++;

print "Creating results table\n";
$status{$n}=system "perl ".$scriptp."createCAGEtable.pl ".$read1.'.CAGEmap '.$gffp.$gff2.' > '.$read1.'.CAGEtable';
die "Problem stage $n" if $status{$n}; $n++;
}


 open (OUT, ">", "CAGE_".$read1.'.r') or die 'could not open CAGE R file';
 print OUT 'source("'.$scriptp.'format_CAGE_table.r")'."\n";
 my $filename=$read1;
 $read1 =~ /\.5mis./;
 $filename=$`;
 print "$filename\n";
 print OUT 'test1=analyseCAGEtable(raw.dir="'.$fastqp.'",cage.dir="./",data="'.$filename.'")'."\n";
 print OUT "q()\n";
 close OUT;


print "formating CAGE data in R.\n";
$status{$n}=system 'R --vanilla --slave < CAGE_'.$read1.'.r > CAGE_'.$read1.'.out';
die "Problem stage $n" if $status{$n}; $n++;

#calculate time taken
my $time1 = new Benchmark;
my $timdiff = timediff($time1, $time0);
print "Running pipeline took ", timestr($timdiff), "\n";
#print "status outputs ", Dumper \%status, "\n\n";











