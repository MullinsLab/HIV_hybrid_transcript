#!/usr/bin/perl -w

##########################################################################################
# Retrieve paired reads from fastq file of original reads
# Author: Wenjie Deng
# Date: 2021-03-04
########################################################################################## 

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my $usage = "usage: retrievePairReadsFromFastq.pl R1Fastq R2Fastq outputPairedR2Fastq\n";

my $r1file = shift or die $usage;
my $r2file = shift or die $usage;
my $outfile = shift or die $usage;

my (%namestatus);
my $r1reads = my $count = my $r2reads = 0;

open R1, $r1file or die "couldn't open $r1file: $!\n";
while (my $name = <R1>) {
	$name =~ s/\R$//;
	my $partialname = $name;
	$partialname =~ /^(\S+)/;
	$partialname = $1;
	my $seq = <R1>;
	my $plus = <R1>;
	my $qual = <R1>;
	if (!$namestatus{$partialname}) {
		$namestatus{$partialname} = 1;
		++$r1reads;
	}else {
		die "duplicate name: $partialname\n";
	}
}
close R1;

open R2, $r2file or die "couldn't open $r2file: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $name = <R2>) {
	$name =~ s/\R$//;
	my $partialname = $name;
	$partialname =~ /^(\S+)/;
	$partialname = $1;
	my $seq = <R2>;
	$seq =~ s/\R$//;
	my $plus = <R2>;
	$plus =~ s/\R$//;
	my $qual = <R2>;
	$qual =~ s/\R$//;
	if ($namestatus{$partialname}) {
		print OUT "$name\n$seq\n$plus\n$qual\n";
		++$r2reads;
	}
	++$count;
}
close R2;
close OUT;

print "retrievePairReadsFromFastq.pl: total $r1reads in $r1file, $count reads in $r2file, paired $r2reads reads to $outfile.\n";
