#!/usr/bin/perl -w

use strict;

my $usage = "perl parseSam.pl inputSamFile outFile mapt2humanLengthCutoff(default: 30)\n";
my $inSam = shift or die $usage;
my $outFile = shift or die $usage;
my $lCut = shift || 30;
my $totalCount = my $validCount = my $invalidCount = my $pairCount = my $mappedCount = my $unmapCount = 0;
my $fwdCount = my $rvsCount = my $gtCutCount = my $lsCutCount = my $fwdGtCutCount = my $rvsGtCutCount = 0;
my $fwdLsCutCount = my $rvsLsCutCount = my $totalAlignNtCount = 0;
my %pairStatus = my %invalidFlag = my %posCoverage = ();
my $refLen = 0;

my $is = my $bp = my $mappingflag = 0;
my $r1ref = my $r2ref = my $r1name = my $r2name = my $orientation = "";
my %chrIsOrientation = my %nameChrIsOrientation = my %flagstatus = ();
open SAM, $inSam or die "couldn't open $inSam: $!\n";

while (my $line = <SAM>) {
	chomp $line;
	next if ($line =~ /^\s*$/ || $line =~ /^@/); 
	++$totalCount;
	my @fields = ();
	@fields = split /\t/, $line;
	my $name = $fields[0];
	my $flag = $fields[1];
	my $ref = $fields[2];
	my $startPos = $fields[3];
	my $pattern = $fields[5];
	my $dir = $fields[8];
	my $seq = $fields[9];
	my $ntCount = my $alignNtCount = my $refAlignNtCount = 0;
	
	
	if ($flag == 0) {	# R1 forward orientation
#		if ($pattern =~ /^\d+M/) {
			my $mapping = countmapping($pattern);
#			if ($mapping >= 30) {
				$is = $startPos;
				$r1ref = $ref;
				$r1name = $name;
				$mappingflag = 1;
				$orientation = "+";
#			}
#		}
	}elsif ($flag == 16) { # R1 reverse complement
#		if ($pattern =~ /\d+M$/) {
			my $mapping = countmapping($pattern);
#			if ($mapping >= 30) {
				$is = $startPos + $mapping - 1;
				$r1ref = $ref;
				$r1name = $name;
				$mappingflag = 1;
				$orientation = "-";
#			}
#		}
	}else {
		if (!$flagstatus{$flag}) {
			++$flagstatus{$flag};
			print "flag: $flag\n";
		}
	}
	
	if ($mappingflag) {		
#		if (!$chrIsOrientation{$r1ref}{$is}{$orientation}) {
			$nameChrIsOrientation{$r1ref}{$is}{$orientation}{$r1name} = $pattern;
#				print OUT "$r1name,$r1ref,$is,$orientation,\n";
# 		}
		++$chrIsOrientation{$r1ref}{$is}{$orientation};		
		
		$r1ref = $r2ref = $r1name = $r2name = $orientation = "";
		$is = $bp = $mappingflag = 0;
	}
		
}
close SAM;

open OUT, ">", $outFile or die "couldn't open $outFile: $!\n";
print OUT "read,reference,IS,orientation,pattern,count\n";
foreach my $ref (sort {$a cmp $b} keys %nameChrIsOrientation) {
	foreach my $is (sort {$a <=> $b} keys %{$nameChrIsOrientation{$ref}}) {		
		foreach my $orientation (keys %{$nameChrIsOrientation{$ref}{$is}}) {
			foreach my $name (keys %{$nameChrIsOrientation{$ref}{$is}{$orientation}}) {
				print OUT "$name,$ref,$is,$orientation,$nameChrIsOrientation{$ref}{$is}{$orientation}{$name},$chrIsOrientation{$ref}{$is}{$orientation}\n";
			}
		}
	}
}
close OUT;


sub countmapping {
	my $pat = shift;
	my $mappingcount = 0;
	my @nums = split /[A-Z]/, $pat;
	my @pats = split /\d+/, $pat;
#	print scalar @nums,', ',scalar @pats,"\n";
	die "impossible: @nums, @pats\n" if (scalar @nums != scalar @pats - 1);
#	print @nums, "\n";
#	print @pats, "\n";
	for (my $i = 0; $i < @pats; $i++) {
		unless ($pats[$i] eq '') {
			if ($pats[$i] eq 'M') {
				$mappingcount += $nums[$i-1];
			}elsif ($pats[$i] eq 'D') {
				$mappingcount += $nums[$i-1];
			}
		}			
	}
	return $mappingcount;
}