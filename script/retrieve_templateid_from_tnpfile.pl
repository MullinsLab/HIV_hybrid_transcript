#!/usr/bin/perl

use strict;

my $usage = "perl retrieve_templateid_from_tnpfile.pl inputtnpfile outfile\n";
my $infile = shift or die $usage;
my $outfile = shift or die $usage;
my %templateName = my %templateCount = ();
my $name = '';
my $totalcount = my $gt2count = my $seqcount = 0;
open IN, $infile or die "coultn't open $infile: $!\n";
while (my $line = <IN>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^>(\S+)/) {
		$name = $1;
	}else {
		push @{$templateName{$line}}, $name;
		++$templateCount{$line};
	}
}
close IN;

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
print OUT "templateid\tcount\treads\n";
foreach my $templateid (sort{$templateCount{$b} <=> $templateCount{$a}} keys %templateCount) {
	print OUT "$templateid\t$templateCount{$templateid}\t", join(",", @{$templateName{$templateid}}), "\n";
	if ($templateCount{$templateid} >= 3) {
		++$gt2count;
	}
	++$totalcount;
	$seqcount += scalar @{$templateName{$templateid}};
}
close OUT;

print "\n* retrieve_templateid_from_tnpfile.pl: total $seqcount reads, $totalcount template ids, $gt2count templates' family size >= 3\n";
