#!/usr/bin/perl -w

##########################################################################################
# Extracts R1 reads starting with the 2nd round primer 
# within first 10bp of the reads and trim the primer, then extract the paired reads,
# trim the 2nd round primer in R2 reads, output both primers trimmed R1 and R2 reads
# Author: Wenjie Deng
# Date: 2017-09-25
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: trimR1PrimerLTR.pl R2Fastq R2Templates outputFile\n";

my $r2file = shift or die $usage;
my $r2templatefile = shift or die $usage;
my $outfile = shift or die $usage;

my (%namestatus, %trimlinkernameSeq, %trimlinkernameQual, %trimalltidnames, %tidcount);
my $trimr1primercount = my $trimLTRcount = my $trimLTRlesscount = my $trimLTRfarcount = 0;
my $lesscount = my $farcount = my $count = my $tidreadcount = my $validtidcount = 0;
open TEMP, $r2templatefile or die "couldn't open $r2templatefile: $!\n";
while (my $line = <TEMP>) {
	chomp $line;
	next if ($line =~ /^\s*$/ or $line =~ /^tid/);
	my @fields = split /\t/, $line;
	my $tid = $fields[0];
	my $counts = $fields[1];
	my $namestring = $fields[2];
#	if (length $tid == 12 and $counts >= 3) {
	if ($tid eq "CGCCCGCTTCTG") {
		my @names = split /,/, $namestring;
		foreach my $name (@names) {
			if (!$namestatus{$name}) {
				$namestatus{$name} = $tid;
				++$validtidcount;
			}else {
				die "duplicate name: $name\n";
			}			
		}
		last;
	}
}
close TEMP;

open R2, $r2file or die "couldn't open $r2file: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $name = <R2>) {
	chomp $name;
	$name =~ /^(\S+)/;
	$name = $1;
	my $seq = <R2>;
	chomp $seq;
	my $plus = <R2>;
	chomp $plus;
	my $qual = <R2>;
	chomp $qual;
	if ($namestatus{$name}) {
		print OUT ">$name\n$seq\n";
		++$count;
	}
}
close R2;
close OUT;

print "$count reads retrieved\n";

