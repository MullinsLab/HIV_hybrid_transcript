#!/usr/bin/perl -w

use strict;

my $usage = "perl parseSamISBreakpoint.pl inputSamFile outFile ltr(3 or 5) mapt2humanLengthCutoff(default: 30)\n";
my $inSam = shift or die $usage;
my $outFile = shift or die $usage;
my $ltr = shift or die $usage;
my $lCut = shift or die $usage;
my $totalCount = my $is = my $bp = my $r1mappingflag = my $r2mappingflag = my $r1mappinglen = my $r2mappinglen = my $r1len = my $r2len = 0;
my $r1ref = my $r2ref = my $r1name = my $r2name = my $orientation = my $r1md = my $r2md = my $r1seq = my $r2seq = my $r1pattern = my $r2pattern = "";
my (%chrIsBpOrientation, %nameChrIsBpOrientation, %flagstatus, %r1nameChrIsBpOrientationidentity, %r2nameChrIsBpOrientationidentity);
my (%r1nameSeq, %r1nameSeqlen, %r1namePattern, %r2nameSeq, %r2nameSeqlen, %r2namePattern);
open SAM, "<", $inSam or die "couldn't open $inSam: $!\n";
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
	my $flen = abs($fields[8]);
	my $seq = $fields[9];
	my $md = $fields[12];
	if ($flag == 99) {	# R1 forward orientation
		$r1mappingflag = $r2mappingflag = $is = $r1len = $r1mappinglen = 0;
		$r1ref = $r1name = $orientation = $r1md = $r1seq = $r1pattern = "";
		if ($pattern =~ /^\d+M/) {
			my $mapping = countmapping($pattern);
			if ($mapping >= 30) {
				$is = $startPos;
				$r1ref = $ref;
				$r1name = $name;
				$r1len = $flen;
				$r1md = $md;
				$r1mappinglen = $mapping;
				$r1seq = $seq;
				$r1pattern = $pattern;
				$r1mappingflag = 1;
				if ($ltr == 3) {
					$orientation = "+";
				}elsif ($ltr == 5) {
					$orientation = "-";
				}else {
					die "Not correct LTR: $ltr\n";
				}				
			}
		}
	}elsif ($flag == 83) { # R1 reverse complement
		$r1mappingflag = $r2mappingflag = $is = $r1len = $r1mappinglen = 0;
		$r1ref = $r1name = $orientation = $r1md = $r1seq = $r1pattern = "";
		if ($pattern =~ /\d+M$/) {
			my $mapping = countmapping($pattern);
			if ($mapping >= 30) {
				$is = $startPos + $mapping - 1;
				$r1ref = $ref;
				$r1name = $name;
				$r1len = $flen;
				$r1md = $md;
				$r1mappinglen = $mapping;
				$r1seq = $seq;
				$r1pattern = $pattern;
				$r1mappingflag = 1;
				if ($ltr == 3) {
					$orientation = "-";
				}elsif ($ltr == 5) {
					$orientation = "+";
				}else {
					die "Not correct LTR: $ltr\n";
				}
			}
		}
	}elsif ($flag == 147) { # R2 reverse complement
		$r2ref = $r2name = $r2md = $r2seq = $r2pattern = "";
		$bp = $r2len = $r2mappingflag = 0;
		if ($pattern =~ /\d+M$/) {
			my $mapping = countmapping($pattern);
			if ($mapping >= 30) {
				$bp = $startPos + $mapping - 1;
				$r2ref = $ref;
				$r2name = $name;
				$r2len = $flen;
				$r2md = $md;
				$r2mappinglen = $mapping;
				$r2seq = $seq;
				$r2pattern = $pattern;
				$r2mappingflag = 1;
			}
		}
	}elsif ($flag == 163) { # R2 forward
		$r2ref = $r2name = $r2md = $r2seq = $r2pattern = "";
		$bp = $r2len = $r2mappingflag = 0;
		if ($pattern =~ /^\d+M/) {
			my $mapping = countmapping($pattern);
			if ($mapping >= 30) {
				$bp = $startPos;
				$r2ref = $ref;
				$r2name = $name;
				$r2len = $flen;
				$r2md = $md;
				$r2mappinglen = $mapping;
				$r2seq = $seq;
				$r2pattern = $pattern;
				$r2mappingflag = 1;
			}
		}
	}else {
		if (!$flagstatus{$flag}) {
			++$flagstatus{$flag};
			print "flag: $flag\n";
		}
	}

	if ($r1mappingflag and ($flag == 147 or $flag == 163)) {
#	if ($r1mappingflag and $r2mappingflag) {
		if ($r1ref eq $r2ref and $r1name eq $r2name and $r1len == $r2len and $r1len <= 1000) {
			$nameChrIsBpOrientation{$r1ref}{$is}{$bp}{$orientation}{$r1name} = 1;
			++$chrIsBpOrientation{$r1ref}{$is}{$bp}{$orientation};
			$r1nameSeq{$r1name} = $r1seq;
			$r2nameSeq{$r2name} = $r2seq;
			$r1nameSeqlen{$r1name} = length $r1seq;
			$r2nameSeqlen{$r2name} = length $r2seq;
			$r1namePattern{$r1name} = $r1pattern;
			$r2namePattern{$r2name} = $r2pattern;
			$r1md =~ s/MD:Z://;
			$r2md =~ s/MD:Z://;	
			my @r1matches = split /[A-Z\^]/, $r1md;
			my @r2matches = split /[A-Z\^]/, $r2md;
			my $r1matchlen = my $r2matchlen = 0;
			foreach my $match (@r1matches) {
				unless ($match eq '') {
					if ($match =~ /^\d+$/) {
						$r1matchlen += $match;
					}else {
						die "MD not formmatted: $line\n";
					}
				}				
			}
			foreach my $match (@r2matches) {
				unless ($match eq '') {
					if ($match =~ /^\d+$/) {
						$r2matchlen += $match;
					}else {
						die "MD not formmatted: $line\n";
					}
				}
			}
			$r1nameChrIsBpOrientationidentity{$r1ref}{$is}{$bp}{$orientation}{$r1name} = int ($r1matchlen / $r1mappinglen * 1000 + 0.5) / 1000;
			$r2nameChrIsBpOrientationidentity{$r1ref}{$is}{$bp}{$orientation}{$r1name} = int ($r2matchlen / $r2mappinglen * 1000 + 0.5) / 1000;
		}
		$r1ref = $r2ref = $r1name = $r2name = $orientation = $r1md = $r2md = $r1seq = $r2seq = $r1pattern = $r2pattern = "";
		$is = $bp = $r1mappingflag = $r2mappingflag = $r1len = $r2len = $r1mappinglen = $r2mappinglen = 0;
	}	
}
close SAM;

open OUT, ">", $outFile or die "couldn't open $outFile: $!\n";
print OUT "read,reference,IS,breakpoint,orientation,r1identity,r2identity,count,r1seq,r1len,r1pattern,r2seq,r2len,r2pattern\n";
foreach my $ref (sort {$a cmp $b} keys %nameChrIsBpOrientation) {
	foreach my $is (sort {$a <=> $b} keys %{$nameChrIsBpOrientation{$ref}}) {
		foreach my $bp (sort {$a <=> $b} keys %{$nameChrIsBpOrientation{$ref}{$is}}) {
			foreach my $orientation (keys %{$nameChrIsBpOrientation{$ref}{$is}{$bp}}) {
				foreach my $name (keys %{$nameChrIsBpOrientation{$ref}{$is}{$bp}{$orientation}}) {					
					print OUT "$name,$ref,$is,$bp,$orientation,$r1nameChrIsBpOrientationidentity{$ref}{$is}{$bp}{$orientation}{$name},$r2nameChrIsBpOrientationidentity{$ref}{$is}{$bp}{$orientation}{$name},$chrIsBpOrientation{$ref}{$is}{$bp}{$orientation},";
					print OUT "$r1nameSeq{$name},$r1nameSeqlen{$name},$r1namePattern{$name},$r2nameSeq{$name},$r2nameSeqlen{$name},$r2namePattern{$name}\n";
				}
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
	die "impossible here: @nums, @pats\n" if (scalar @nums != scalar @pats - 1);
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