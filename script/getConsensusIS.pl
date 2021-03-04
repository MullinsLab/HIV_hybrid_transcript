#!/usr/bin/perl -w

##########################################################################################
# From IS list by parsing .sam file and R1_tid.txt file, get the consensus IS for each  
# template
# Author: Wenjie Deng
# Date: 2017-12-07
# Modified: 2017-12-25
# added gene location, orientation, start and end positions for each IS
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: getConsensusIS.pl ISList R1tidReadnamefile humanGeneGffFile outputTemplateConsensusIS\n";

my $isfile = shift or die $usage;
my $tidfile = shift or die $usage;
my $gfffile = shift or die $usage;
my $outfile = shift or die $usage;

my (%namestatus, %nameIS, %nameRef, %nameDir, %chromoGene);
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
	my $dir = $fields[3];
	$nameIS{$name} = $is;
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
print OUT "templateID,is,chr,chr_orientation,gene,gene_orientation,gene_start,gene_end,R1_count,R1_map2human,total_R1_reads\n";
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
				++$isCount{$ref}{$dir}{$is};
			}			
		}
		if ($actualcount >= 3) {
			my $maxRef = my $maxDir = "";
			my $maxIS = my $maxcount = 0;
			foreach my $ref (keys %isCount) {
				foreach my $dir (keys %{$isCount{$ref}}) {
					foreach my $is (keys %{$isCount{$ref}{$dir}}) {
						my $count = $isCount{$ref}{$dir}{$is};
						if ($count > $maxcount) {
							$maxcount = $count;
							$maxIS = $is;
							$maxRef = $ref;
							$maxDir = $dir;
						}						
					}
				}
			}
			my $fraction = $maxcount / $actualcount;
			if ($fraction > 0.5) {
				my $isgene = my $dir = "N/A";
				my $genestart = my $geneend = "N/A";
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
				print OUT "$tid,$maxIS,$maxRef,$maxDir,$isgene,$dir,$genestart,$geneend,$maxcount,$actualcount,$counts\n";
			}
		}		
	}	
}
close TID;
close OUT;

