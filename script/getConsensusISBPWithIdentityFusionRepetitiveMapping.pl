#!/usr/bin/perl -w

##########################################################################################
# From IS list with mapping identities by parsing .sam file, get the 
# consensus IS and breackpoint passing the cutoff of identity (default 0.99)
# Author: Wenjie Deng
# Date: 2021-05-05
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: getConsensusISBPWithIdentityFusionRepetitiveMapping.pl ISList humanGeneGffFile outputConsensusISbreakpointCSVFile identityCutoff\n";

my $isfile = shift or die $usage;
my $gfffile = shift or die $usage;
my $outfile = shift or die $usage;
my $cutoff = shift || 0.99;
my $fusionfile = my $repetitivefile = $outfile;
$fusionfile =~ s/\.csv/_fusion.csv/;
$repetitivefile =~ s/\.csv/_repetitive.csv/;

my (%namestatus, %nameIS, %nameBP, %nameRef, %nameDir, %chromoGene, %refbpdirmulti, %passcutoffrefisbpdirmulti, %fusionpasscutoffrefisbpdirmulti, %repetitivepasscutoffrefisbpdirmulti);
my $lesscount = my $farcount = my $count = my $tidreadcount = 0;
open IS, $isfile or die "couldn't open $isfile: $!\n";
while (my $line = <IS>) {
	chomp $line;
	next if ($line =~ /^\s*$/ or $line =~ /^read,/);
	my @fields = split /,/, $line;
	my $name = $fields[0];
	my $r2ref = $fields[1];
	my $bp = $fields[2]; 
	my $r2flag = $fields[3];
	my $r1ref = $fields[7];
	my $r1pos = $fields[8];
	my $r1flag = $fields[9];
	my $r1mapq = $fields[10];
	my $r1as = $fields[11];
	my $r1xs = $fields[12];	
	my $dir = $fields[13];	
	my $r1identity = $fields[14];
	my $r2identity = $fields[15];
	my $multi = $fields[16];
	my $r1pattern = $fields[19];
	my $r1mlen = $fields[20];
	my $r2pattern = $fields[23];
	
	$refbpdirmulti{$r2ref}{$bp}{$dir} = $multi;
	if (!$namestatus{$name}) {
		$namestatus{$name} = 1;
	}else {
		die "duplicate name: $name\n";
	}
	
	if ($r1identity >= $cutoff and $r2identity >= $cutoff) {
		if ($r1flag == 99 and $r2flag == 147) {	# R1 forward and properly mapped
			if ($r1pattern =~ /^\d+M/ and $r2pattern =~ /\d+M$/) {
				my $is = $r1pos;
				if ($r1mapq > 0) {	# not repetitive						
					++$passcutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}else {
					++$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}
			}
		}elsif ($r1flag == 97 and $r2flag == 145) { # R1 forward and fusion
			if ($r1pattern =~ /^\d+M/ and $r2pattern =~ /\d+M$/) {
				my $is = $r1pos;
				if ($r1mapq > 0) {	# not repetitive						
					++$fusionpasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}else {
					++$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}
			}
		}elsif ($r1flag == 83 and $r2flag == 163) {	# R1 reverse and properly mapped
			if ($r1pattern =~ /\d+M$/ and $r2pattern =~ /^\d+M/) {
				my $is = $r1pos + $r1mlen - 1;
				if ($r1mapq > 0) {	# not repetitive						
					++$passcutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}else {
					++$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}
			}
		}elsif ($r1flag == 81 and $r2flag == 161) {	# R1 reverse and fusion
			if ($r1pattern =~ /\d+M$/ and $r2pattern =~ /^\d+M/) {
				my $is = $r1pos + $r1mlen - 1;
				if ($r1mapq > 0) {	# not repetitive						
					++$fusionpasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}else {
					++$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}
			}
		}		
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
print OUT "ISchr,IS,BPchr,BP,Chr_orientation,Gene,Gene_orientation,Gene_start,Gene_end,total_BP_count,pass_IS_BP_identity_".$cutoff."_count\n";
foreach my $r1ref (sort {$a cmp $b} keys %passcutoffrefisbpdirmulti) {
	foreach my $r2ref (sort {$a cmp $b} keys %{$passcutoffrefisbpdirmulti{$r1ref}}) {
		foreach my $is (sort {$a <=> $b} keys %{$passcutoffrefisbpdirmulti{$r1ref}{$r2ref}}) {
			foreach my $bp (sort {$a <=> $b} keys %{$passcutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}}) {
				foreach my $dir (keys %{$passcutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}}) {
					my $isgene = my $genedir = "NA";
					my $genestart = my $geneend = "NA";
					foreach my $gene (keys %{$chromoGene{$r1ref}}) {
						if ($is >= $chromoGene{$r1ref}{$gene}{start} and $is <= $chromoGene{$r1ref}{$gene}{end}) {
							$isgene = $gene;
							$genestart = $chromoGene{$r1ref}{$gene}{start};
							$geneend = $chromoGene{$r1ref}{$gene}{end};
							if ($chromoGene{$r1ref}{$gene}{dir} eq "+") {
								$genedir = $dir;
							}elsif ($chromoGene{$r1ref}{$gene}{dir} eq "-") {
								if ($dir eq "+") {
									$genedir = "-";
								}else {
									$genedir = "+";
								}
							}else {
								die "No gene orientation for $r1ref, $gene\n";
							}
							last;
						}
					}
					print OUT "$r1ref,$is,$r2ref,$bp,$dir,$isgene,$genedir,$genestart,$geneend,$refbpdirmulti{$r2ref}{$bp}{$dir},$passcutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir}\n";
				}
			}
		}
	}
	
}
close OUT;

open FUSION, ">", $fusionfile or die "couldn't open $fusionfile: $!\n";
print FUSION "ISchr,IS,BPchr,BP,Chr_orientation,Gene,Gene_orientation,Gene_start,Gene_end,total_BP_count,pass_IS_BP_identity_".$cutoff."_count\n";
foreach my $r1ref (sort {$a cmp $b} keys %fusionpasscutoffrefisbpdirmulti) {
	foreach my $r2ref (sort {$a cmp $b} keys %{$fusionpasscutoffrefisbpdirmulti{$r1ref}}) {
		foreach my $is (sort {$a <=> $b} keys %{$fusionpasscutoffrefisbpdirmulti{$r1ref}{$r2ref}}) {
			foreach my $bp (sort {$a <=> $b} keys %{$fusionpasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}}) {
				foreach my $dir (keys %{$fusionpasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}}) {
					my $isgene = my $genedir = "NA";
					my $genestart = my $geneend = "NA";
					foreach my $gene (keys %{$chromoGene{$r1ref}}) {
						if ($is >= $chromoGene{$r1ref}{$gene}{start} and $is <= $chromoGene{$r1ref}{$gene}{end}) {
							$isgene = $gene;
							$genestart = $chromoGene{$r1ref}{$gene}{start};
							$geneend = $chromoGene{$r1ref}{$gene}{end};
							if ($chromoGene{$r1ref}{$gene}{dir} eq "+") {
								$genedir = $dir;
							}elsif ($chromoGene{$r1ref}{$gene}{dir} eq "-") {
								if ($dir eq "+") {
									$genedir = "-";
								}else {
									$genedir = "+";
								}
							}else {
								die "No gene orientation for $r1ref, $gene\n";
							}
							last;
						}
					}
					print FUSION "$r1ref,$is,$r2ref,$bp,$dir,$isgene,$genedir,$genestart,$geneend,$refbpdirmulti{$r2ref}{$bp}{$dir},$fusionpasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir}\n";
				}
			}
		}
	}
	
}
close FUSION;

open REP, ">", $repetitivefile or die "couldn't open $repetitivefile: $!\n";
print REP "ISchr,IS,BPchr,BP,Chr_orientation,Gene,Gene_orientation,Gene_start,Gene_end,total_BP_count,pass_IS_BP_identity_".$cutoff."_count\n";
foreach my $r1ref (sort {$a cmp $b} keys %repetitivepasscutoffrefisbpdirmulti) {
	foreach my $r2ref (sort {$a cmp $b} keys %{$repetitivepasscutoffrefisbpdirmulti{$r1ref}}) {
		foreach my $is (sort {$a <=> $b} keys %{$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}}) {
			foreach my $bp (sort {$a <=> $b} keys %{$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}}) {
				foreach my $dir (keys %{$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}}) {
					my $isgene = my $genedir = "NA";
					my $genestart = my $geneend = "NA";
					foreach my $gene (keys %{$chromoGene{$r1ref}}) {
						if ($is >= $chromoGene{$r1ref}{$gene}{start} and $is <= $chromoGene{$r1ref}{$gene}{end}) {
							$isgene = $gene;
							$genestart = $chromoGene{$r1ref}{$gene}{start};
							$geneend = $chromoGene{$r1ref}{$gene}{end};
							if ($chromoGene{$r1ref}{$gene}{dir} eq "+") {
								$genedir = $dir;
							}elsif ($chromoGene{$r1ref}{$gene}{dir} eq "-") {
								if ($dir eq "+") {
									$genedir = "-";
								}else {
									$genedir = "+";
								}
							}else {
								die "No gene orientation for $r1ref, $gene\n";
							}
							last;
						}
					}
					print REP "$r1ref,$is,$r2ref,$bp,$dir,$isgene,$genedir,$genestart,$geneend,$refbpdirmulti{$r2ref}{$bp}{$dir},$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir}\n";
				}
			}
		}
	}
	
}
close REP;
