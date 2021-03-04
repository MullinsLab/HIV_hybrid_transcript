#!/usr/bin/perl -w

##########################################################################################
# Retrieve IS and breakpoint associated R1 and R2 (only human part) reads 
# Author: Wenjie Deng
# Date: 2019-06-24
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: retrieveISR1R2Reads.pl R1Fastq R2Fastq bwaParsedFile outputR1Fastq outputR2Fastq\n";

my $r1file = shift or die $usage;
my $r2file = shift or die $usage;
my $csvfile = shift or die $usage;
my $outr1file = shift or die $usage;
my $outr2file = shift or die $usage;

my (%namestatus);
my $r1reads = my $count = my $r2reads = 0;
open CSV, $csvfile or die "couldn't open $csvfile: $!\n";
while (my $line = <CSV>) {
	chomp $line;
	next if ($line =~ /^\s*$/ or $line =~ /^read/);
	my @fields = split /,/, $line;
	my $name = "@".$fields[0];
	$namestatus{$name} = 1;
	++$count;
}
close CSV;

open R1, $r1file or die "couldn't open $r1file: $!\n";
open OUT, ">", $outr1file or die "couldn't open $outr1file: $!\n";
while (my $name = <R1>) {
	chomp $name;
	my $partialname = $name;
	$partialname =~ /^(\S+)/;
	$partialname = $1;
	my $seq = <R1>;
	chomp $seq;
	my $plus = <R1>;
	chomp $plus;
	my $qual = <R1>;
	chomp $qual;
	if ($namestatus{$partialname}) {
		print OUT "$name\n$seq\n$plus\n$qual\n";
		++$r1reads;
	}
}
close R1;
close OUT;

open R2, $r2file or die "couldn't open $r2file: $!\n";
open OUT, ">", $outr2file or die "couldn't open $outr2file: $!\n";
while (my $name = <R2>) {
	chomp $name;
	my $partialname = $name;
	$partialname =~ /^(\S+)/;
	$partialname = $1;
	my $seq = <R2>;
	chomp $seq;
	my $plus = <R2>;
	chomp $plus;
	my $qual = <R2>;
	chomp $qual;
	if ($namestatus{$partialname}) {
		print OUT "$name\n$seq\n$plus\n$qual\n";
		++$r2reads;
	}
}
close R2;
close OUT;

print "Total $count reads in $csvfile file. Retrieved $r1reads R1 reads, $r2reads R2 reads.\n";

