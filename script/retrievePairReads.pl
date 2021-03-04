#!/usr/bin/perl -w

##########################################################################################
# Extracts paired reads from read fastq file 
# Author: Wenjie Deng
# Date: 2019-05-30
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: retrievePairReads.pl R1Fastq R2Templates outputR1Fastq\n";

my $r1file = shift or die $usage;
my $r2templatefile = shift or die $usage;
my $outfile = shift or die $usage;

my (%namestatus);
my $r1reads = my $count = my $validtidcount = my $validtidreadcount = 0;
open TEMP, $r2templatefile or die "couldn't open $r2templatefile: $!\n";
while (my $line = <TEMP>) {
	chomp $line;
	next if ($line =~ /^\s*$/ or $line =~ /^templateid/);
	my @fields = split /\t/, $line;
	my $tid = $fields[0];
	my $counts = $fields[1];
	my $namestring = $fields[2];
	if (length $tid == 12 and $counts >= 3) {
		++$validtidcount;
		my @names = split /,/, $namestring;
		foreach my $name (@names) {
			if (!$namestatus{$name}) {
				$namestatus{$name} = $tid;
				++$validtidreadcount;
			}else {
				die "duplicate name: $name\n";
			}			
		}
	}
}
close TEMP;

open R1, $r1file or die "couldn't open $r1file: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
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
	++$count;
}
close R1;
close OUT;

print "\nvalid template Id (12nt, count >= 3) in $r2templatefile: $validtidcount, $validtidreadcount reads\n";
print "total $count reads in $r1file, paired $r1reads reads to $outfile.\n";
