#!/usr/bin/perl

##########################################################################################
# program: template_primer_fasta.pl
# From before and after trimmed files to locate trimmed sequences, output as fasta file
# usage: perl template_primer_fasta.pl before_trim_fastq after_trim_fastq
# Author: Wenjie Deng
# Date: 2016-05-04
##########################################################################################


use strict;

my $usage = "usage: perl template_primer_fasta.pl before_trim_fastq after_trim_fastq 3or5(3' or 5')\n";
my $beforefile = shift or die $usage;
my $afterfile = shift or die $usage;
my $flank = shift || 5;
my $outfile = $afterfile;
if ($flank == 3) {
	$outfile =~ s/\.fastq$/_3tnp\.fasta/;
}elsif ($flank == 5) {
	$outfile =~ s/\.fastq$/_5tnp\.fasta/;
}else {
	die "couldn't recognize: $flank\n";
}

my $count = my $fastacount = 0;
my %afternameSeq = ();
open AFTER, $afterfile or die "couldn't open $afterfile: $!\n";
while (my $name = <AFTER>) {
	chomp $name;
	my $seq = <AFTER>;
	chomp $seq;
	$afternameSeq{$name} = $seq;
	my $plus = <AFTER>;
	my $qual = <AFTER>;
	++$count;
}
close AFTER;

open BEFORE, $beforefile or die "couldn't open $beforefile: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $name = <BEFORE>) {
	chomp $name;
	my $seq = <BEFORE>;
	chomp $seq;
	my $plus = <BEFORE>;
	my $qual = <BEFORE>;
	if ($afternameSeq{$name}) {
		my $afterseq = $afternameSeq{$name};
		my $idx = index($seq, $afterseq);
		my $tnpseq = "";
		if ($flank == 5) {
			$tnpseq = substr($seq, 0, $idx);		
		}else {
			my $afterlen = length $afterseq;
			$tnpseq = substr($seq, $idx+$afterlen);
		}		
		print OUT ">$name\n$tnpseq\n";
		++$fastacount;
	}
}
close BEFORE;
close OUT;
print "\n* template_primer_fasta.pl: $flank' total $count reads after, write $fastacount to fasta\n\n";