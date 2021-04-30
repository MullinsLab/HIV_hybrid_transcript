#!/usr/bin/perl

use strict;
use v5.10;

my $ltrfile = shift;
my $adptfile = shift;
my $linkerfile = shift;
my $ltr = my $adpt = my $linker = "";
open LTR, $ltrfile;
while (my $line = <LTR>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$ltr = $line;
		last;
	}
}
close LTR;

open ADPT, $adptfile;
while (my $line = <ADPT>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$adpt = $line;
		last;
	}
}
close ADPT;

open LINKER, $linkerfile;
while (my $line = <LINKER>) {
	$line =~ s/\R$//;
	next if $line =~ /^\s*$/;
	unless ($line =~ /^>/) {
		$linker = $line;
		last;
	}
}
close LINKER;

my $ltrlen = length $ltr;
my $adptlen = length $adpt;
my $linkerlen = length $linker;

my $ltrminoverlap = int($ltrlen * 0.9 + 0.5);
my $adptminoverlap = int($adptlen * 0.9 + 0.5);
my $linkerminoverlap = int($linkerlen * 0.9 + 0.5);

print "$ltr\n$adpt\n$linker\n$ltrminoverlap\n$adptminoverlap\n$linkerminoverlap\n";
