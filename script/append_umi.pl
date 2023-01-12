#!/usr/bin/perl

# from _UMI.fasta file, append reads' umi into _R2_xxx_trimmed.csv

use strict;
use File::Copy;

my $usage = "perl append_umi.pl inputUMIfile inputR2TrimmedCSVfile outputR2TrimmedPlusUMIcsvfile\n";
my $umifile = shift or die $usage;
my $csvfile = shift or die $usage;
my $outfile = shift or die $usage;
my %nameUmi = ();
my $name = '';
my $totalcount = my $gt2count = my $seqcount = 0;
open IN, $umifile or die "coultn't open $umifile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
	}else {
		$nameUmi{$name} = $line;
	}
}
close IN;

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
open CSV, "<", $csvfile or die "couldn't open $csvfile: $!\n";
while (my $line = <CSV>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	my ($name, $trimmed, $seq) = split /,/, $line;
	if ($nameUmi{$name}) {
		print OUT "$line,$nameUmi{$name}\n";
	}else {
		print "*** No UMI for R2 read $name ***\n";
		print OUT "$line,NA\n";
	}
}
close CSV;
close OUT;

