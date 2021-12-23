#!/usr/bin/perl -w

##########################################################################################
# From IS list with mapping identities by parsing .sam file, get the 
# consensus IS and breackpoint passing the cutoff of identity (default 0.99)
# Author: Wenjie Deng
# Date: 2021-05-05
# Modified: calculate R1 human consensus for each IS and BP combination, the IS is the
# beginning of sequence
# Date: 2021-11-19
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: getConsensusISBPWithIdentityFusionRepetitiveMapping.pl ISList humanGeneGffFile outputConsensusISbreakpointCSVFile identityCutoff\n";

my $isfile = shift or die $usage;
my $gfffile = shift or die $usage;
my $outfile = shift or die $usage;
my $cutoff = shift || 0.99;
my $fusionfile = my $repetitivefile = my $collapsedfile = $outfile;
$fusionfile =~ s/\.csv/_fusion.csv/;
$repetitivefile =~ s/\.csv/_repetitive.csv/;

my (%namestatus, %nameIS, %nameBP, %nameRef, %nameDir, %chromoGene, %refbpdirmulti, %refisbpdirmulti, %passcutoffrefisbpdirmulti, %fusionpasscutoffrefisbpdirmulti, %repetitivepasscutoffrefisbpdirmulti, %refisbpdirr1humanseqs, %refisbpdirr2humanseqs);
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
	my $r1seq = $fields[18];
	my $r1pattern = $fields[20];
	my $r1mlen = $fields[21];
	my $r2seq = $fields[23];
	my $r2pattern = $fields[25];
	my $r2mlen = $fields[26];
	my $is = 0;
	$refbpdirmulti{$r2ref}{$bp}{$dir} = $multi;
	if (!$namestatus{$name}) {
		$namestatus{$name} = 1;
	}else {
		die "duplicate name: $name\n";
	}
	
	if ($r1flag == 99 and $r2flag == 147) {	# R1 forward and properly mapped
		if ($r1pattern =~ /^\d+M/ and $r2pattern =~ /\d+M$/) {
			$is = $r1pos;
			++$refisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
			if ($r1identity >= $cutoff and $r2identity >= $cutoff) {
				if ($r1mapq > 0) {	# not repetitive						
					++$passcutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}else {
					++$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}
				my $r1humanseq = substr($r1seq, 0, $r1mlen);
				push @{$refisbpdirr1humanseqs{$r1ref}{$r2ref}{$is}{$bp}{$dir}}, $r1humanseq;
				my $r2humanseq = substr($r2seq, 0, $r2mlen);
				push @{$refisbpdirr2humanseqs{$r1ref}{$r2ref}{$is}{$bp}{$dir}}, $r2humanseq;	
			}			
		}
	}elsif ($r1flag == 97 and $r2flag == 145) { # R1 forward and fusion
		if ($r1pattern =~ /^\d+M/ and $r2pattern =~ /\d+M$/) {
			$is = $r1pos;
			++$refisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
			if ($r1identity >= $cutoff and $r2identity >= $cutoff) {
				if ($r1mapq > 0) {	# not repetitive						
					++$fusionpasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}else {
					++$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}
				my $r1humanseq = substr($r1seq, 0, $r1mlen);
				push @{$refisbpdirr1humanseqs{$r1ref}{$r2ref}{$is}{$bp}{$dir}}, $r1humanseq;
				my $r2humanseq = substr($r2seq, 0, $r2mlen);
				push @{$refisbpdirr2humanseqs{$r1ref}{$r2ref}{$is}{$bp}{$dir}}, $r2humanseq;
			}			
		}
	}elsif ($r1flag == 83 and $r2flag == 163) {	# R1 reverse and properly mapped
		if ($r1pattern =~ /\d+M$/ and $r2pattern =~ /^\d+M/) {
			$is = $r1pos + $r1mlen - 1;
			++$refisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
			if ($r1identity >= $cutoff and $r2identity >= $cutoff) {
				if ($r1mapq > 0) {	# not repetitive						
					++$passcutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}else {
					++$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}
				my $r1humanseq = substr($r1seq, 0, $r1mlen);
				push @{$refisbpdirr1humanseqs{$r1ref}{$r2ref}{$is}{$bp}{$dir}}, $r1humanseq;
				my $r2humanseq = substr($r2seq, 0, $r2mlen);
				push @{$refisbpdirr2humanseqs{$r1ref}{$r2ref}{$is}{$bp}{$dir}}, $r2humanseq;
			}			
		}
	}elsif ($r1flag == 81 and $r2flag == 161) {	# R1 reverse and fusion
		if ($r1pattern =~ /\d+M$/ and $r2pattern =~ /^\d+M/) {
			$is = $r1pos + $r1mlen - 1;
			++$refisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
			if ($r1identity >= $cutoff and $r2identity >= $cutoff) {
				if ($r1mapq > 0) {	# not repetitive						
					++$fusionpasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}else {
					++$repetitivepasscutoffrefisbpdirmulti{$r1ref}{$r2ref}{$is}{$bp}{$dir};
				}
				my $r1humanseq = substr($r1seq, 0, $r1mlen);
				push @{$refisbpdirr1humanseqs{$r1ref}{$r2ref}{$is}{$bp}{$dir}}, $r1humanseq;
				my $r2humanseq = substr($r2seq, 0, $r2mlen);
				push @{$refisbpdirr2humanseqs{$r1ref}{$r2ref}{$is}{$bp}{$dir}}, $r2humanseq;
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

if (%passcutoffrefisbpdirmulti) {
	outputConsensusISbreakpointCSVFile($outfile, \%passcutoffrefisbpdirmulti, \%refisbpdirr1humanseqs, \%refisbpdirr2humanseqs, \%refisbpdirmulti, \%refbpdirmulti, \%chromoGene);
}
if (%fusionpasscutoffrefisbpdirmulti) {
	outputConsensusISbreakpointCSVFile($fusionfile, \%fusionpasscutoffrefisbpdirmulti, \%refisbpdirr1humanseqs, \%refisbpdirr2humanseqs, \%refisbpdirmulti, \%refbpdirmulti, \%chromoGene);
}
if (%repetitivepasscutoffrefisbpdirmulti) {
	outputConsensusISbreakpointCSVFile($repetitivefile, \%repetitivepasscutoffrefisbpdirmulti, \%refisbpdirr1humanseqs, \%refisbpdirr2humanseqs, \%refisbpdirmulti, \%refbpdirmulti, \%chromoGene);
}


sub outputConsensusISbreakpointCSVFile {
	my ($outfile, $passcutoffrefisbpdirmulti_ref, $refisbpdirr1humanseqs_ref, $refisbpdirr2humanseqs_ref, $refisbpdirmulti_ref, $refbpdirmulti_ref, $chromoGene_ref) = @_;
	open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
	print OUT "ISchr,IS,BPchr,BP,Chr_orientation,Gene,Gene_orientation,Gene_start,Gene_end,total_BP_count,total_IS_BP_count,pass_IS_BP_identity_".$cutoff."_count,IS_consensus,BP_consensus\n";
	foreach my $r1ref (sort {$a cmp $b} keys %{$passcutoffrefisbpdirmulti_ref}) {
		foreach my $r2ref (sort {$a cmp $b} keys %{$passcutoffrefisbpdirmulti_ref->{$r1ref}}) {
			foreach my $is (sort {$a <=> $b} keys %{$passcutoffrefisbpdirmulti_ref->{$r1ref}->{$r2ref}}) {
				foreach my $bp (sort {$a <=> $b} keys %{$passcutoffrefisbpdirmulti_ref->{$r1ref}->{$r2ref}->{$is}}) {
					foreach my $dir (keys %{$passcutoffrefisbpdirmulti_ref->{$r1ref}->{$r2ref}->{$is}->{$bp}}) {
						my $isgene = my $genedir = "NA";
						my $genestart = my $geneend = "NA";
						foreach my $gene (keys %{$chromoGene_ref->{$r1ref}}) {
							if ($is >= $chromoGene_ref->{$r1ref}->{$gene}->{start} and $is <= $chromoGene_ref->{$r1ref}->{$gene}->{end}) {
								$isgene = $gene;
								$genestart = $chromoGene_ref->{$r1ref}->{$gene}->{start};
								$geneend = $chromoGene_ref->{$r1ref}->{$gene}->{end};
								if ($chromoGene_ref->{$r1ref}->{$gene}->{dir} eq "+") {
									$genedir = $dir;
								}elsif ($chromoGene_ref->{$r1ref}->{$gene}->{dir} eq "-") {
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
						my $iscollapsedfile = $is."_".$bp."_iscollapsed.fasta";
						my $isconsensus = consensus_seqs($iscollapsedfile, \@{$refisbpdirr1humanseqs_ref->{$r1ref}->{$r2ref}->{$is}->{$bp}->{$dir}});
						$isconsensus =~ s/-//g;
						my $bpcollapsedfile = $is."_".$bp."_bpcollapsed.fasta";
						my $bpconsensus = consensus_seqs($bpcollapsedfile, \@{$refisbpdirr2humanseqs_ref->{$r1ref}->{$r2ref}->{$is}->{$bp}->{$dir}});
						$bpconsensus =~ s/-//g;
						print OUT "$r1ref,$is,$r2ref,$bp,$dir,$isgene,$genedir,$genestart,$geneend,$refbpdirmulti_ref->{$r2ref}->{$bp}->{$dir},$refisbpdirmulti_ref->{$r1ref}->{$r2ref}->{$is}->{$bp}->{$dir},$passcutoffrefisbpdirmulti_ref->{$r1ref}->{$r2ref}->{$is}->{$bp}->{$dir},$isconsensus,$bpconsensus\n";
					}
				}
			}
		}	
	}
	close OUT;
}


sub get_consensus {
	my $nt = shift;
	my %consnt = (
		'-' => '-',
		'A' => 'A',
		'C' => 'C',
		'G' => 'G',
		'T' => 'T',
		'AG' => 'R',
		'CT' => 'Y',
		'CG' => 'S',
		'AT' => 'W',
		'GT' => 'K',
		'AC' => 'M',
		'CGT' => 'B',
		'ACG' => 'V',
		'AGT' => 'D',
		'ACT' => 'H',
		'ACGT' => 'N',
	);
	return $consnt{$nt};
}


sub consensus_seqs {
	my $file = shift;
	my $seqs_ref = shift;
	my %seqCount = ();
	my $consensus = "";
	foreach my $seq (@{$seqs_ref}) {
		if (!$seqCount{$seq}) {
			$seqCount{$seq} = 0;
		}
		++$seqCount{$seq};
	}
	if (keys %seqCount == 1) {
		($consensus) = keys %seqCount;
		my ($value) = values %seqCount;
	}else {
		my @sortedseqs = sort {$seqCount{$b} <=> $seqCount{$a}} keys %seqCount;
		my $idx = my $count = 0;
		foreach my $seq (@sortedseqs) {
			$count += $seqCount{$seq};
		}
		if ($seqCount{$sortedseqs[0]} > 0.5 * $count) {
			$consensus = $sortedseqs[0];
		}else {
			open TMP, ">", $file or die "couldn't open $file: $!\n";
			foreach my $seq (@sortedseqs) {
				++$idx;
				my $name = $idx."_".$seqCount{$seq};
				print TMP ">$name\n$seq\n";
			}
			my $align_file = $file;
			$align_file =~ s/\.fasta/_align.fasta/;
			system("muscle -quiet -in $file -out $align_file");
			my $name = "";
			my $alignlen = 0;
			my (@names, %nameSeq, %nameStart, %nameEnd, %nameNts, %posCount, %posNtCount);
			open ALIGN, "<", $align_file or die "couldn't open $align_file: $!\n";
			while (my $line = <ALIGN>) {
				chomp $line;
				if ($line =~ /^>(\S+)/) {
					$name = $1;
					push @names, $name;
				}else {
					$nameSeq{$name} .= $line;
				}
			}
			close ALIGN;
			unlink $file;
			unlink $align_file;
			foreach my $name (@names) {
				my $seq = $nameSeq{$name};
				my @nts = split //, $seq;
				@{$nameNts{$name}} = split //, $seq;
				if (!$alignlen) {
					$alignlen = length $seq;
				}
				if (length $seq != $alignlen) {
					die "sequences are not aligned\n";
				}
				for (my $i = 0; $i < $alignlen; $i++) {
					if ($nts[$i] =~ /[ACGT]/) {
						$nameStart{$name} = $i;
						last;
					}
				}
				for (my $i = $alignlen - 1; $i >= 0; $i--) {
					if ($nts[$i] =~ /[ACGT]/) {
						$nameEnd{$name} = $i;
						last;
					}
				}
			}
			for (my $i = 0; $i < $alignlen; $i++) {
				foreach my $name (@names) {
					my $duplicates = 0;
					if ($name =~ /_(\d+)$/) {
						$duplicates = $1;
					}else {
						die "name not formatted: $name\n";
					}
					if ($i >= $nameStart{$name} and $i <= $nameEnd{$name}) {
						$posCount{$i} += $duplicates;
						$posNtCount{$i}{$nameNts{$name}[$i]} += $duplicates;
					}
				}
				if ($posCount{$i} >= 0.5 * $count) {
					# get consensus Nt
					my $consnt = '';
					my @sortedNts = sort{$posNtCount{$i}{$b} <=> $posNtCount{$i}{$a}} keys %{$posNtCount{$i}};
					if (scalar @sortedNts == 1 or $posNtCount{$i}{$sortedNts[0]} > $posNtCount{$i}{$sortedNts[1]}) {
						$consnt = $sortedNts[0];
					}elsif ($posNtCount{$i}{$sortedNts[0]} == $posNtCount{$i}{$sortedNts[1]}) {						
						if ($sortedNts[0] eq '-') {
							$consnt = $sortedNts[1];
						}else {
							$consnt = $sortedNts[0].$sortedNts[1];
						}
						if ($sortedNts[2] and $posNtCount{$i}{$sortedNts[2]} == $posNtCount{$i}{$sortedNts[0]}) {
							$consnt .= $sortedNts[2];
						}
						if ($sortedNts[3] and $posNtCount{$i}{$sortedNts[3]} == $posNtCount{$i}{$sortedNts[0]}) {
							$consnt .= $sortedNts[3];
						}
						if ($sortedNts[4] and $posNtCount{$i}{$sortedNts[4]} == $posNtCount{$i}{$sortedNts[0]}) {
							$consnt .= $sortedNts[4];
						}
					}else {
						die "impossible! $posNtCount{$i}{$sortedNts[0]} < $posNtCount{$i}{$sortedNts[1]}\n";
					}
					if (length $consnt > 1) {
						my @nts = split //, $consnt;
						my @sortnts = sort {$a cmp $b} @nts;
						$consnt = join('', @sortnts);
					}
					$consensus .= get_consensus($consnt)
				}else {
					last;
				}
			} 
		}
	}		
	return $consensus;
}
