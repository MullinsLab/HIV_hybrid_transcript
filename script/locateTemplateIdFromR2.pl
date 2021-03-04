#!/usr/bin/perl -w

##########################################################################################
# Extracts R2 reads starting with the 2nd round primer 
# within first 10bp of the reads and trim the primer, then trim "GTTATGGTACTT" to locate
# template Id, output template Id associated with read names 
# Author: Wenjie Deng
# Date: 2017-09-27
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: locateTemplateIdFromR2.pl R2_fastq outfile\n";

my $r2File = shift or die $usage;
my $outfile = shift or die $usage;
my $r2Primer = "CCGCTCCGTCCGACGACTCACTATA";
my $linker = "GTTATGGTACTT";
my (%r2trimprimernameSeq, %r2tidnames, %r2tidcount);
my $probename = '';
my $count = my $r2trimprimercount = my $r2linkercount = my $linkerfarcount = my $farcount = my $r2lesscount = my $r2farcount = 0;

open R2, $r2File or die "couldn't open $r2File: $!\n";
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
	my $idx = index($seq, $r2Primer);
	if ($idx >= 0) { # there is a hit
		if ($idx < 10) {
			my $trimStart = $idx + length($r2Primer); 
			my $trimseq = substr($seq, $trimStart);
			if ($trimseq) {
				$r2trimprimernameSeq{$name} = $trimseq;
				++$r2trimprimercount;
			}
		}else {
			++$farcount;
		}
	}
	++$count;
}
close R2;

foreach my $name (keys %r2trimprimernameSeq) {
	my $seq = $r2trimprimernameSeq{$name};
	my $idx = index($seq, $linker);
	if ($idx >= 0) { # there is a hit
		if ($idx < 20) {
			my $trimStart = $idx + length($linker); 
			my $tid = substr($seq, 0, $idx);
			if ($tid) {
				push @{$r2tidnames{$tid}}, $name;
				++$r2tidcount{$tid};
				++$r2linkercount;
			}
		}else {
			++$linkerfarcount;
		}
	}
}

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
print OUT "tid\tcounts\tnames\n";
foreach my $tid (sort {$r2tidcount{$b} <=> $r2tidcount{$a}} keys %r2tidcount) {
	print OUT "$tid\t$r2tidcount{$tid}\t";
	foreach my $name (@{$r2tidnames{$tid}}) {
		print OUT "$name,";
	}
	print OUT "\n";
}

print "total $count R2 reads, trim 2nd primer $r2trimprimercount, farcount: $farcount, trim linker $r2linkercount, farcount: $linkerfarcount\n";

