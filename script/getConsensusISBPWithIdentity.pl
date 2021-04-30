#!/usr/bin/perl -w

##########################################################################################
# From IS list with mapping identities by parsing .sam file, get the 
# consensus IS and breackpoint passing the cutoff of identity (default 0.99)
# Author: Wenjie Deng
# Date: 2021-03-05
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: getConsensusISBPWithIdentity.pl ISList humanGeneGffFile outputConsensusISbreakpointCSVFile identityCutoff\n";

my $isfile = shift or die $usage;
my $gfffile = shift or die $usage;
my $outfile = shift or die $usage;
my $cutoff = shift || 0.99;

my (%namestatus, %nameIS, %nameBP, %nameRef, %nameDir, %chromoGene, %tidinfo, %refisbpdirmulti, %passcutoffrefisbpdirmulti);
my $lesscount = my $farcount = my $count = my $tidreadcount = 0;
open IS, $isfile or die "couldn't open $isfile: $!\n";
while (my $line = <IS>) {
	chomp $line;
	next if ($line =~ /^\s*$/ or $line =~ /read,reference/);
	my @fields = split /,/, $line;
	my $name = $fields[0];
	my $ref = $fields[1];
	my $is = $fields[2]; 
	my $bp = $fields[3];
	my $dir = $fields[4];
	my $r1identity = $fields[5];
	my $r2identity = $fields[6];
	my $multi = $fields[7];
	$nameIS{$name} = $is;
	$nameBP{$name} = $bp;
	$nameRef{$name} = $ref;
	$nameDir{$name} = $dir;
	$refisbpdirmulti{$ref}{$is}{$bp}{$dir} = $multi;
	$namestatus{$name} = 1;
	if ($r1identity >= $cutoff and $r2identity >= $cutoff) {
		++$passcutoffrefisbpdirmulti{$ref}{$is}{$bp}{$dir};
	}
}
close IS;

open GFF, $gfffile or die "couldn't open $gfffile: $!\n";
while (my $line = <GFF>) {
	chomp $line;
	next if $line =~ /^#/;
	my @fields = split /\t/, $line;
	if ($fields[2] eq "gene") {
		my $gene = '';
		my $geneId = 0;
		my $attribute = $fields[8];
		if ($attribute =~ /;Name=(.*?);/) {
			$gene = $1;
		}
		if ($attribute =~ /GeneID:(\d+)/) {
			$geneId = $1;
		}
		$chromoGene{$fields[0]}{$gene}{start} = $fields[3];
		$chromoGene{$fields[0]}{$gene}{end} = $fields[4];
		$chromoGene{$fields[0]}{$gene}{dir} = $fields[6];
		$chromoGene{$fields[0]}{$gene}{id} = $geneId;
	}
}
close GFF;

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
print OUT "Chr,IS,BP,Chr_orientation,Gene,Gene_orientation,Gene_start,Gene_end,total_count,pass_identity_".$cutoff."_count\n";
foreach my $ref (sort {$a cmp $b} keys %passcutoffrefisbpdirmulti) {
	foreach my $is (sort {$a <=> $b} keys %{$passcutoffrefisbpdirmulti{$ref}}) {
		foreach my $bp (sort {$a <=> $b} keys %{$passcutoffrefisbpdirmulti{$ref}{$is}}) {
			foreach my $dir (keys %{$passcutoffrefisbpdirmulti{$ref}{$is}{$bp}}) {
				my $isgene = my $genedir = "NA";
				my $genestart = my $geneend = "NA";
				foreach my $gene (keys %{$chromoGene{$ref}}) {
					if ($is >= $chromoGene{$ref}{$gene}{start} and $is <= $chromoGene{$ref}{$gene}{end}) {
						$isgene = $gene;
						$genestart = $chromoGene{$ref}{$gene}{start};
						$geneend = $chromoGene{$ref}{$gene}{end};
						if ($chromoGene{$ref}{$gene}{dir} eq "+") {
							$genedir = $dir;
						}elsif ($chromoGene{$ref}{$gene}{dir} eq "-") {
							if ($dir eq "+") {
								$genedir = "-";
							}else {
								$genedir = "+";
							}
						}else {
							die "No gene orientation for $ref, $gene\n";
						}
						last;
					}
				}
				print OUT "$ref,$is,$bp,$dir,$isgene,$genedir,$genestart,$geneend,$refisbpdirmulti{$ref}{$is}{$bp}{$dir},$passcutoffrefisbpdirmulti{$ref}{$is}{$bp}{$dir}\n";
			}
		}
	}
}
close OUT;
