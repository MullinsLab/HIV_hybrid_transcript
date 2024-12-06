#!/usr/bin/perl -w

use strict;
use warnings;
use v5.10;

my $usage = "perl parseSamBreakpoint.pl inputSamFile outFile R1trimmecsvfile R2trimmecsvfile ltr(3 or 5) map2humanLengthCutoff(default: 20)\n";
my $inSam = shift or die $usage;
my $outFile = shift or die $usage;
my $r1trfile = shift or die $usage;
my $r2trfile = shift or die $usage;
my $ltr = shift or die $usage;
my $lCut = shift or die $usage;
my $totalCount = my $is = my $bp = my $r1mappingflag = my $r2mappingflag = my $r1mappinglen = my $r2mappinglen = my $r1len = my $r2len = my $r1pos = my $r1flag = my $r2flag = my $r1as = my $r2as = my $r1xs = my $r2xs = 0;
my $r1ref = my $r2ref = my $r1name = my $r2name = my $orientation = my $r1md = my $r2md = my $r1seq = my $r2seq = my $r1pattern = my $r2pattern = my $r1mapq = my $r2mapq = "";
my (%chrIsBpOrientation, %nameChrIsBpOrientation, %flagstatus, %r1nameChrIsBpOrientationidentity, %r2nameChrIsBpOrientationidentity, %r1nameAs, %r2nameAs, %r1nameXs, %r2nameXs);
my (%r1nameSeq, %r1nameSeqlen, %r1namePattern, %r2nameSeq, %r2nameSeqlen, %r2namePattern, %r1nameRef, %r1namePos, %r1nameFlag, %r2nameFlag, %r1nameMapq, %r2nameMapq, %r1nameMaplen, %r2nameMaplen);
my (%r1nameTrimmed, %r1nameRemain, %r2nameTrimmed, %r2nameRemain, %r2nameUmi);
open R1, "<", $r1trfile or die "couldn't open $r1trfile: $!\n";
while (my $line = <R1>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/);
	my ($name, $trimmed, $remain) = split /,/, $line;
	$r1nameTrimmed{$name} = $trimmed;
	$r1nameRemain{$name} = $remain;
}
close R1;
open R2, "<", $r2trfile or die "couldn't open $r2trfile: $!\n";
while (my $line = <R2>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/);
	my ($name, $trimmed, $remain, $umi) = split /,/, $line;
	$r2nameTrimmed{$name} = $trimmed;
	$r2nameRemain{$name} = $remain;
	$r2nameUmi{$name} = $umi;
}
close R2;

open SAM, "<", $inSam or die "couldn't open $inSam: $!\n";
while (my $line = <SAM>) {
	$line =~ s/\R$//;
	next if ($line =~ /^\s*$/ || $line =~ /^@/); 
	++$totalCount;
	my @fields = ();
	@fields = split /\t/, $line;
	my $name = $fields[0];
	my $flag = $fields[1];
	my $ref = $fields[2];
	my $startPos = $fields[3];
	my $mapq = $fields[4];
	my $pattern = $fields[5];
	my $flen = abs($fields[8]);
	my $seq = $fields[9];
	my $md = $fields[12];
	my $as = $fields[14];
	my $xs = $fields[15];
	if ($flag == 99 or $flag == 97) {	# R1 forward orientation
		$r1mappingflag = $r2mappingflag = $r1pos = $r1len = $r1mappinglen = $r1flag = $r1as = $r1xs = 0;
		$r1ref = $r1name = $orientation = $r1md = $r1seq = $r1pattern = $r1mapq = "";
		my $mapping = countmapping($pattern);
		if ($mapping >= $lCut) {
			$r1pos = $startPos;
			$r1ref = $ref;
			$r1name = $name;
			$r1flag = $flag;
			$r1mapq = $mapq;
			$r1len = $flen;
			$r1md = $md;
			$r1as = $as;
			$r1as =~ s/AS\:i\://;
			$r1xs = $xs;
			$r1xs =~ s/XS\:i\://;
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
	}elsif ($flag == 83 or $flag == 81) { # R1 reverse complement
		$r1mappingflag = $r2mappingflag = $r1pos = $r1len = $r1mappinglen = $r1flag = $r1as = $r1xs = 0;
		$r1ref = $r1name = $orientation = $r1md = $r1seq = $r1pattern = $r1mapq = "";
		my $mapping = countmapping($pattern);
		if ($mapping >= $lCut) {
			$r1pos = $startPos;
			$r1ref = $ref;
			$r1name = $name;
			$r1flag = $flag;
			$r1mapq = $mapq;
			$r1len = $flen;
			$r1md = $md;
			$r1as = $as;
			$r1as =~ s/AS\:i\://;
			$r1xs = $xs;
			$r1xs =~ s/XS\:i\://;
			$r1mappinglen = $mapping;
			$r1seq = reversecomplement($seq);
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
	}elsif ($flag == 147 or $flag == 145) { # R2 reverse complement
		$r2ref = $r2name = $r2md = $r2seq = $r2pattern = $r2mapq = "";
		$bp = $r2len = $r2mappingflag = $r2flag = $r2as = $r2xs = 0;
		if ($pattern =~ /\d+M$/) {
			my $mapping = countmapping($pattern);
			if ($mapping >= $lCut) {
				$bp = $startPos + $mapping - 1;
				$r2ref = $ref;
				$r2name = $name;
				$r2flag = $flag;
				$r2mapq = $mapq;
				$r2len = $flen;
				$r2md = $md;
				$r2as = $as;
				$r2as =~ s/AS\:i\://;
				$r2xs = $xs;
				$r2xs =~ s/XS\:i\://;
				$r2mappinglen = $mapping;
				$r2seq = reversecomplement($seq);
				$r2pattern = $pattern;
				$r2mappingflag = 1;
			}
		}
	}elsif ($flag == 163 or $flag == 161) { # R2 forward
		$r2ref = $r2name = $r2md = $r2seq = $r2pattern = $r2mapq = "";
		$bp = $r2len = $r2mappingflag = $r2flag = $r2as = $r2xs = 0;
		if ($pattern =~ /^\d+M/) {
			my $mapping = countmapping($pattern);
			if ($mapping >= $lCut) {
				$bp = $startPos;
				$r2ref = $ref;
				$r2name = $name;
				$r2flag = $flag;
				$r2mapq = $mapq;
				$r2len = $flen;
				$r2md = $md;
				$r2as = $as;
				$r2as =~ s/AS\:i\://;
				$r2xs = $xs;
				$r2xs =~ s/XS\:i\://;
				$r2mappinglen = $mapping;
				$r2seq = $seq;
				$r2pattern = $pattern;
				$r2mappingflag = 1;
			}
		}
	}else {
		if (!$flagstatus{$flag}) {
			++$flagstatus{$flag};
			#print "flag: $flag\n";
		}
	}

	if ($r1mappingflag and ($flag == 147 or $flag == 163 or $flag == 145 or $flag == 161)) {
		if ($r1name eq $r2name) {
			$nameChrIsBpOrientation{$r2ref}{$bp}{$orientation}{$r2name} = 1;
			++$chrIsBpOrientation{$r2ref}{$bp}{$orientation};
			if ($r1seq ne $r1nameRemain{$r1name}) {
				die "R1 mapping to human are different: $r1name: $1seq vs $r1nameRemain{$r1name}\n";
			}
			if ($r2seq ne $r2nameRemain{$r2name}) {
				die "R2 mapping to human are different: $r2name: $2seq vs $r2nameRemain{$r2name}\n";
			}
			$r1nameSeq{$r1name} = $r1seq;
			$r2nameSeq{$r2name} = $r2seq;
			$r1nameSeqlen{$r1name} = length $r1seq;
			$r2nameSeqlen{$r2name} = length $r2seq;
			$r1namePattern{$r1name} = $r1pattern;
			$r2namePattern{$r2name} = $r2pattern;
			$r1nameRef{$r1name} = $r1ref;
			$r1namePos{$r1name} = $r1pos;
			$r1nameFlag{$r1name} = $r1flag;
			$r2nameFlag{$r2name} = $r2flag;
			$r1nameMapq{$r1name} = $r1mapq;
			$r2nameMapq{$r2name} = $r2mapq;
			$r1nameMaplen{$r1name} = $r1mappinglen;
			$r2nameMaplen{$r2name} = $r2mappinglen;
			$r1nameAs{$r1name} = $r1as;
			$r2nameAs{$r2name} = $r2as;
			$r1nameXs{$r1name} = $r1xs;
			$r2nameXs{$r2name} = $r2xs;
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
			$r1nameChrIsBpOrientationidentity{$r1ref}{$bp}{$orientation}{$r1name} = int ($r1matchlen / $r1mappinglen * 1000 + 0.5) / 1000;
			$r2nameChrIsBpOrientationidentity{$r2ref}{$bp}{$orientation}{$r2name} = int ($r2matchlen / $r2mappinglen * 1000 + 0.5) / 1000;
		}
		$r1ref = $r2ref = $r1name = $r2name = $orientation = $r1md = $r2md = $r1seq = $r2seq = $r1pattern = $r2pattern = $r1mapq = $r2mapq = "";
		$r1pos = $bp = $r1mappingflag = $r2mappingflag = $r1len = $r2len = $r1mappinglen = $r2mappinglen = $r1flag = $r2flag = $r1as = $r2as = $r1xs = $r2xs = 0;
	}	
}
close SAM;

open OUT, ">", $outFile or die "couldn't open $outFile: $!\n";
print OUT "read,r2ref,breakpoint,r2flag,r2mapq,r2as,r2xs,r1ref,r1pos,r1flag,r1mapq,r1as,r1xs,orientation,r1identity,r2identity,count,r1ltr,r1seq,r1len,r1pattern,r1mlen,r2linker,r2seq,r2len,r2pattern,r2mlen,umi\n";
foreach my $ref (sort {$a cmp $b} keys %nameChrIsBpOrientation) {
	foreach my $bp (sort {$a <=> $b} keys %{$nameChrIsBpOrientation{$ref}}) {
		foreach my $orientation (keys %{$nameChrIsBpOrientation{$ref}{$bp}}) {
			foreach my $name (keys %{$nameChrIsBpOrientation{$ref}{$bp}{$orientation}}) {					
				print OUT "$name,$ref,$bp,$r2nameFlag{$name},$r2nameMapq{$name},$r2nameAs{$name},$r2nameXs{$name},$r1nameRef{$name},$r1namePos{$name},$r1nameFlag{$name},$r1nameMapq{$name},";
				print OUT "$r1nameAs{$name},$r1nameXs{$name},$orientation,$r1nameChrIsBpOrientationidentity{$r1nameRef{$name}}{$bp}{$orientation}{$name},$r2nameChrIsBpOrientationidentity{$ref}{$bp}{$orientation}{$name},";
				print OUT "$chrIsBpOrientation{$ref}{$bp}{$orientation},$r1nameTrimmed{$name},$r1nameSeq{$name},$r1nameSeqlen{$name},$r1namePattern{$name},$r1nameMaplen{$name},$r2nameTrimmed{$name},$r2nameSeq{$name},$r2nameSeqlen{$name},$r2namePattern{$name},$r2nameMaplen{$name},$r2nameUmi{$name}\n";
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

sub reversecomplement {
	my $seq = shift;
	$seq = reverse $seq;
	$seq =~ tr /ACGTacgt/TGCAtgca/;	
	return $seq;
}