#!/usr/bin/perl -w

##########################################################################################
# Retrieve sequence of trimmed parts after cutting the adapter and saved as .csv file
# Author: Wenjie Deng
# Date: 2021-05-07
########################################################################################## 

use strict;
use warnings;
use v5.10;
use Getopt::Long;
use File::Basename;

my $usage = "usage: retrieveTrimmedSeq.pl beforeTrimmingFastq afterTrimmingFastq outputTrimmingFasta\n";

my $btfile = shift or die $usage;
my $atfile = shift or die $usage;
my $outfile = shift or die $usage;

my (%namestatus, %btnameseq, %atnameseq);
my $btreads = my $outcount = my $atreads = 0;

open BT, $btfile or die "couldn't open $btfile: $!\n";
while (my $name = <BT>) {
	$name =~ s/\R$//;
	$name =~ /^(\S+)/;
	$name = $1;
	$name =~ s/^\@//g;
	my $seq = <BT>;
	$seq =~ s/\R$//;
	my $plus = <BT>;
	my $qual = <BT>;
	if (!$namestatus{$name}) {
		$namestatus{$name} = 1;
		++$btreads;
	}else {
		die "duplicate name: $name\n";
	}
	$btnameseq{$name} = $seq;
}
close BT;

open AT, $atfile or die "couldn't open $atfile: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
while (my $name = <AT>) {
	$name =~ s/\R$//;
	$name =~ /^(\S+)/;
	$name = $1;
	$name =~ s/^\@//g;
	my $seq = <AT>;
	$seq =~ s/\R$//;
	my $plus = <AT>;
	my $qual = <AT>;
	if ($namestatus{$name}) {
		if ($seq) {
			my $idx = index($btnameseq{$name}, $seq);
			if ($idx > 0) {
				my $trimmedseq = substr($btnameseq{$name}, 0, $idx);
				print OUT "$name,$trimmedseq,$seq\n";
				++$outcount;
			}elsif ($idx == 0) {
				die "read $name not trimmed: $seq vs $btnameseq{$name}\n";
			}else {
				die "No matches of $seq to $btnameseq{$name}\n"; 
			}
		}else { # no sequence left, all trimmed
			print OUT "$name,$btnameseq{$name}\n";
			++$outcount;
		}		
	}else {
		die "No name $name in before trimmed fastq file\n";
	}
	++$atreads;
}
close AT;
close OUT;

print "retrieveTrimmedSeq.pl: total $btreads before trimming, $atreads after trimming, $outcount trimmed sequences to $outfile.\n";
