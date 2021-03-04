#!/usr/bin/perl -w

##########################################################################################
# From Trimmed 3'LTR R1 fastq file and R2 templates file, make R1 templates and fastq files
# Author: Wenjie Deng
# Date: 2019-06-04
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: makeR1Templates.pl R1Fastq R2Templates outputR1Templates outR1TemplateFastq\n";

my $r1file = shift or die $usage;
my $r2templatefile = shift or die $usage;
my $outfile = shift or die $usage;
my $outr1file = shift or die $usage;

my (%namestatus, %r1templatenamestatus);
my $r1reads = my $count = my $validtidcount = my $r1tidcount = 0;
open R1, $r1file or die "couldn't open $r1file: $!\n";
while (my $name = <R1>) {
	chomp $name;
	my $partialname = $name;
	$partialname =~ /^(\S+)/;
	$partialname = $1;
	my $seq = <R1>;
	my $plus = <R1>;
	my $qual = <R1>;
	if (!$namestatus{$partialname}) {
		$namestatus{$partialname} = 1;
		++$r1reads;
	}
}
close R1;

open TEMP, $r2templatefile or die "couldn't open $r2templatefile: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $line = <TEMP>) {
	chomp $line;
	next if ($line =~ /^\s*$/ or $line =~ /^templateid/);
	my @fields = split /\t/, $line;
	my $tid = $fields[0];
	my $counts = $fields[1];
	my $namestring = $fields[2];
	my @validnames = ();
	my $count = 0;
	if (length $tid == 12) {
		my @names = split /,/, $namestring;
		foreach my $name (@names) {
			if ($namestatus{$name}) {
				push @validnames, $name;
				++$count;
			}		
		}
		++$validtidcount;
	}
	if ($count >= 3) {
		foreach my $name (@validnames) {
			$r1templatenamestatus{$name} = 1;
		}
		print OUT "$tid\t$count\t", join(',', @validnames), "\n";
		++$r1tidcount;
	}
	
}
close TEMP;
close OUT;

open R1, $r1file or die "couldn't open $r1file: $!\n";
open OUT, ">", $outr1file or die "couldn't open $outr1file: $!\n";
while (my $name = <R1>) {
	chomp $name;
	my $partialname = $name;
	$partialname =~ /^(\S+)/;
	$partialname = $1;
	my $seq = <R1>;
	my $plus = <R1>;
	my $qual = <R1>;
	if ($r1templatenamestatus{$partialname}) {
		print OUT "$name\n$seq$plus$qual";
	}
}
close R1;
close OUT;

print "\nR1 tid with valid template Id (12nt, >= 3 reads): $r1tidcount\n";
