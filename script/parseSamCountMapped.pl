#!/usr/bin/perl -w

use strict;

my $usage = "perl parseSam.pl inputSamFile cutoff\n";
my $inSam = shift or die $usage;
my $fCut = shift || 0.6;
my $R1Mapped = my $R2Mapped = my $totalCount = 0;
my $R1MappedGtCut = my $R1MappedLtCut = my $R2MappedGtCut = my $R2MappedLtCut = 0;
my $validCount = my $invalidCount = my $mapCount = my $unmapCount = my $pairCount = 0;
my (%flagstatus, %pairStatus, %invalidFlag);
open SAM, $inSam or die "couldn't open $inSam: $!\n";
while (my $line = <SAM>) {
	chomp $line;
	next if ($line =~ /^\s*$/ || $line =~ /^@/); 
	++$totalCount;
	my @fields = ();
	@fields = split /\t/, $line;
	my $name = $fields[0];
	my $flag = $fields[1];
	my $startPos = $fields[3];
	my $pattern = $fields[5];
	my $dir = $fields[8];
	my $seq = $fields[9];
	my $ntCount = my $alignNtCount = my $refAlignNtCount = 0;
#	print "$fCut\n$name\n$flag\n$pattern\n$seq\n";
	if ($pattern !~ /\d+H/) {
		++$validCount;
		if (!$pairStatus{$name}) {
			$pairStatus{$name} = 0;
		}
		++$pairStatus{$name};
		if ($pairStatus{$name} == 2) {
			++$pairCount;
		}
		die "more than 1 pair: $name\n" if ($pairStatus{$name} > 2);
		if ($pattern eq '*' && $dir == 0) {
			++$unmapCount;
#			if ($flag != 77 && $flag != 141) {
#				print "$line\n";
#			}
		}else {
			++$mapCount;
			if ($flag == 99 or $flag == 83) {	# R1 
				++$R1Mapped;
			}elsif ($flag == 147 or $flag == 163) { # R2 
				++$R2Mapped;
			}			
			my $len = length $seq;
			my @nums = split /[A-Z]/, $pattern;
			my @pats = split /\d+/, $pattern;
#			print scalar @nums,', ',scalar @pats,"\n";
			die "impossible: $name, @nums, @pats\n" if (scalar @nums != scalar @pats - 1);
#			print @nums, "\n";
#			print @pats, "\n";
			for (my $i = 0; $i < @pats; $i++) {
				unless ($pats[$i] eq '') {
					if ($pats[$i] eq 'S') {
						$ntCount += $nums[$i-1];
					}elsif ($pats[$i] eq 'M') {
						$ntCount += $nums[$i-1];
						$alignNtCount += $nums[$i-1];
						$refAlignNtCount += $nums[$i-1];
					}elsif ($pats[$i] eq 'D') {
						$refAlignNtCount += $nums[$i-1];
					}elsif ($pats[$i] eq 'I') {
						$ntCount += $nums[$i-1];
						$alignNtCount += $nums[$i-1];
					}else {
						die "$name, $pattern, other $pats[$i]";
					}
				}			
			}
			die "length not same, $name, $len, $ntCount\n" if ($ntCount != $len);
			my $fraction = $alignNtCount / $ntCount;
			if ($fraction >= $fCut) {			
				if ($flag == 99 or $flag == 83) {	# R1 
					++$R1MappedGtCut;
				}elsif ($flag == 147 or $flag == 163) { # R2 
					++$R2MappedGtCut;
				}
			}else {
				if ($flag == 99 or $flag == 83) {	# R1 
					++$R1MappedLtCut;
				}elsif ($flag == 147 or $flag == 163) { # R2 
					++$R2MappedLtCut;
				}
			}
		}		
	}else {
		++$invalidCount;
		$invalidFlag{$flag} = 1;
#		print "$line\n";
	}
}
close SAM;

print "Total: $totalCount, valid: $validCount, invalid: $invalidCount, Mapped: $mapCount, Unmapped: $unmapCount.\n";
print "Mapped R1: $R1Mapped, great than cutoff $fCut: $R1MappedGtCut, less than cutoff: $R1MappedLtCut.\n";
print "Mapped R2: $R2Mapped, great than cutoff $fCut: $R2MappedGtCut, less than cutoff: $R2MappedLtCut.\n";