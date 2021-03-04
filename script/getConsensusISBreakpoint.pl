#!/usr/bin/perl -w

##########################################################################################
# From IS list by parsing .sam file and R1_sickle_pair_3LTR_templates.txt file, get the 
# consensus IS and breackpoint for each template
# Author: Wenjie Deng
# Date: 2019-06-11
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: getConsensusISBreakpoint.pl ISList R1tidReadnamefile humanGeneGffFile outputTemplateConsensusISbreakpointCSVFile\n";

my $isfile = shift or die $usage;
my $tidfile = shift or die $usage;
my $gfffile = shift or die $usage;
my $outfile = shift or die $usage;

my (%namestatus, %nameIS, %nameBP, %nameRef, %nameDir, %chromoGene, %tidinfo);
my $trimr1primercount = my $trimLTRcount = my $trimLTRlesscount = my $trimLTRfarcount = 0;
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
	$nameIS{$name} = $is;
	$nameBP{$name} = $bp;
	$nameRef{$name} = $ref;
	$nameDir{$name} = $dir;
	$namestatus{$name} = 1;
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

open TID, $tidfile or die "couldn't open $tidfile: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
print OUT "TemplateID,IS,BP,Chr,Chr_orientation,Gene,Gene_orientation,Gene_start,Gene_end,R1_count,R1_map2human,total_R1_reads\n";
while (my $line = <TID>) {
	chomp $line;
	next if ($line =~ /^\s*$/ or $line =~ /^tid/);
	my @fields = split /\t/, $line;
	my $tid = $fields[0];
	my $counts = $fields[1];
	my $namestring = $fields[2];
	my $actualcount = 0;
	my %isCount = ();
	if (length $tid == 12 and $counts >= 3) {
		my @names = split /,/, $namestring;
		foreach my $name (@names) {
			$name =~ s/^\@//;
			if ($namestatus{$name}) {
				++$actualcount;
				my $ref = $nameRef{$name};
				my $dir = $nameDir{$name};
				my $is = $nameIS{$name};
				my $bp = $nameBP{$name};
				++$isCount{$ref}{$dir}{$is}{$bp};
			}			
		}
		if ($actualcount >= 3) {
			my $maxRef = my $maxDir = "";
			my $maxIS = my $maxBP = my $maxcount = 0;
			foreach my $ref (keys %isCount) {
				foreach my $dir (keys %{$isCount{$ref}}) {
					foreach my $is (keys %{$isCount{$ref}{$dir}}) {
						foreach my $bp (keys %{$isCount{$ref}{$dir}{$is}}) {
							my $count = $isCount{$ref}{$dir}{$is}{$bp};
							if ($count > $maxcount) {
								$maxcount = $count;
								$maxIS = $is;
								$maxBP = $bp;							
								$maxRef = $ref;
								$maxDir = $dir;
							}	
						}						
					}
				}
			}
			my $fraction = $maxcount / $actualcount;
			if ($fraction > 0.5) {
				my $isgene = my $dir = "NA";
				my $genestart = my $geneend = "NA";
				foreach my $gene (keys %{$chromoGene{$maxRef}}) {
					if ($maxIS >= $chromoGene{$maxRef}{$gene}{start} && $maxIS <= $chromoGene{$maxRef}{$gene}{end}) {
						$isgene = $gene;
						$genestart = $chromoGene{$maxRef}{$gene}{start};
						$geneend = $chromoGene{$maxRef}{$gene}{end};
						if ($chromoGene{$maxRef}{$gene}{dir} eq "+") {
							$dir = $maxDir;
						}elsif ($chromoGene{$maxRef}{$gene}{dir} eq "-") {
							if ($maxDir eq "+") {
								$dir = "-";
							}else {
								$dir = "+";
							}
						}else {
							die "No gene orientation for $maxRef, $gene\n";
						}
						last;
					}
				}
				print OUT "$tid,$maxIS,$maxBP,$maxRef,$maxDir,$isgene,$dir,$genestart,$geneend,$maxcount,$actualcount,$counts\n";
				my $tidplus = $tid.'('.$counts.','.$actualcount.','.$maxcount.')';
				push @{$tidinfo{$maxRef}{$maxIS}{$maxBP}{$maxDir}{$isgene}{$dir}{$genestart}{$geneend}}, $tidplus;
			}
		}		
	}	
}
close TID;
close OUT;
my $txtout = $outfile;
$txtout =~ s/csv$/txt/;
open OUT, ">", $txtout or die "couldn't open $txtout: $!\n";
print OUT "Chr\tIS\tBP\tChr_orientation\tGene\tGene_orientation\tGene_start\tGene_end\tTIDs(Reads,MappedToHuman,Consensus)\n";
foreach my $ref (sort {$a cmp $b} keys %tidinfo) {
	foreach my $is (sort {$a <=> $b} keys %{$tidinfo{$ref}}) {
		foreach my $bp (sort {$a <=> $b} keys %{$tidinfo{$ref}{$is}}) {
			foreach my $refdir (keys %{$tidinfo{$ref}{$is}{$bp}}) {
				foreach my $gene (keys %{$tidinfo{$ref}{$is}{$bp}{$refdir}}) {
					foreach my $genedir (keys %{$tidinfo{$ref}{$is}{$bp}{$refdir}{$gene}}) {
						foreach my $start (keys %{$tidinfo{$ref}{$is}{$bp}{$refdir}{$gene}{$genedir}}) {
							foreach my $end (keys %{$tidinfo{$ref}{$is}{$bp}{$refdir}{$gene}{$genedir}{$start}}) {
								print OUT "$ref\t$is\t$bp\t$refdir\t$gene\t$genedir\t$start\t$end\t",join(",", @{$tidinfo{$ref}{$is}{$bp}{$refdir}{$gene}{$genedir}{$start}{$end}}),"\n";
							}
						}
					}
				}
				
			}
		}
	}
}
close OUT;



