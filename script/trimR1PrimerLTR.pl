#!/usr/bin/perl -w

##########################################################################################
# Extracts R1 reads starting with the 2nd round primer 
# within first 10bp of the reads and trim the primer, then extract the paired reads,
# trim the 2nd round primer in R2 reads, output both primers trimmed R1 and R2 reads
# Author: Wenjie Deng
# Date: 2017-09-25
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: trimR1PrimerLTR.pl R1Fastq R2Templates outputR1trimR1primerLTRfastq tidReadnamefile\n";

my $r1file = shift or die $usage;
my $r2templatefile = shift or die $usage;
my $outfile = shift or die $usage;
my $tidfile = shift or die $usage;
my $r1Primer = "GCCCGTCTGTTGTGTGACTCTGGTAACTAGAGAT";
my $linker = "TCTCTAGCA";
#my $linker = "CCCTCAGACC";
my (%namestatus, %trimlinkernameSeq, %trimlinkernameQual, %trimalltidnames, %tidcount);
my $trimr1primercount = my $trimLTRcount = my $trimLTRlesscount = my $trimLTRfarcount = 0;
my $lesscount = my $farcount = my $count = my $tidreadcount = my $validtidcount = 0;
open TEMP, $r2templatefile or die "couldn't open $r2templatefile: $!\n";
while (my $line = <TEMP>) {
	chomp $line;
	next if ($line =~ /^\s*$/ or $line =~ /^tid/);
	my @fields = split /\t/, $line;
	my $tid = $fields[0];
	my $counts = $fields[1];
	my $namestring = $fields[2];
	if (length $tid == 12 and $counts >= 3) {
		my @names = split /,/, $namestring;
		foreach my $name (@names) {
			if (!$namestatus{$name}) {
				$namestatus{$name} = $tid;
				++$validtidcount;
			}else {
				die "duplicate name: $name\n";
			}			
		}
	}
}
close TEMP;

open R1, $r1file or die "couldn't open $r1file: $!\n";
while (my $name = <R1>) {
	chomp $name;
	$name =~ /^(\S+)/;
	$name = $1;
	my $seq = <R1>;
	chomp $seq;
	my $plus = <R1>;
	chomp $plus;
	my $qual = <R1>;
	chomp $qual;
	if ($namestatus{$name}) {
		my $idx = index($seq, $r1Primer);
		if ($idx >= 0) { # there is a hit
			if ($idx < 10) {
				my $trimStart = $idx + length($r1Primer); 
				my $trimseq = substr($seq, $trimStart);
				my $trimQual = substr($qual, $trimStart);
				if ($trimseq and length $trimseq >= 50) {
					++$trimr1primercount;
					
					my $linkeridx = index($trimseq, $linker);
					if ($linkeridx >= 0) { # there is a hit
						if ($linkeridx > 20 and $linkeridx < 40) {
#						if ($idx >= 0 and $idx < 10) {
							my $trimStart = $linkeridx + length($linker); 
							my $trimlinkerseq = substr($trimseq, $trimStart);
							my $trimlinkerqual = substr($trimQual, $trimStart);
							if ($trimlinkerseq and length $trimlinkerseq >= 30) {
								$trimlinkernameSeq{$name} = $trimlinkerseq;
								$trimlinkernameQual{$name} = $trimlinkerqual;
								push @{$trimalltidnames{$namestatus{$name}}}, $name;
								++$tidcount{$namestatus{$name}};
								++$trimLTRcount;
							}else {
								++$trimLTRlesscount;
							}
						}else {
							++$trimLTRfarcount;
						}
					}					
				}else {
					++$lesscount;
				}
			}else {
				++$farcount;
			}
		}
		++$tidreadcount;
	}
	++$count;
}
close R1;

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
open TID, ">", $tidfile or die "couldn't open $tidfile: $!\n";
print TID "tid\tcounts\tnames\n";
foreach my $tid (sort {$tidcount{$b} <=> $tidcount{$a}} keys %tidcount) {
	my $counts = scalar @{$trimalltidnames{$tid}};
	if ($counts >= 3) {
		print TID "$tid\t$counts\t";
		foreach my $name (@{$trimalltidnames{$tid}}) {
			print TID "$name,";
			print OUT "$name\n$trimlinkernameSeq{$name}\n\+\n$trimlinkernameQual{$name}\n";
		}
	}
	print TID "\n";
}
close OUT;
close TID;
print "\nR2 reads with valid template Id (12nt, >= 3 reads): $validtidcount\n";
print "total $count R1 reads, $tidreadcount with valid template Id, $trimr1primercount trimmed R1 primer, lesscount $lesscount, farcount $farcount\n";
print "$trimLTRcount trimmed LTR, trimmed LTR lesscount $trimLTRlesscount, trimmed LTR farcount $trimLTRfarcount\n";

