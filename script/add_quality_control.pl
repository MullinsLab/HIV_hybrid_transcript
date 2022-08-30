#!/usr/bin/perl

# from IS .csv file, append quality control column

use strict;
use File::Copy;

my $usage = "perl add_quality_control.pl inputISCSVfile outputIScsvfile\n";
my $csvfile = shift or die $usage;
my $outfile = shift or die $usage;
my (@isbps, %isbpline, %isbppasscutoff, %chrisstatus, %chrbpstatus, %chrisflag, %chrbpflag);
my $header = '';
my $totalcount = my $gt2count = my $seqcount = 0;

open CSV, "<", $csvfile or die "couldn't open $csvfile: $!\n";
while (my $line = <CSV>) {
	chomp $line;
	next if $line =~ /^\s*$/;
	if ($line =~ /^ISchr,IS,/) {
		$header = $line;
	}else {
		my @fields = split /,/, $line;
		my $ischr = $fields[0];
		my $is = $fields[1];
		my $bpchr = $fields[2];
		my $bp = $fields[3];
		my $passcutoff = $fields[12];
		if (!$chrisstatus{$ischr}{$is}) {
			$chrisstatus{$ischr}{$is} = 0;
		}
		++$chrisstatus{$ischr}{$is};
		if (!$chrbpstatus{$bpchr}{$bp}) {
			$chrbpstatus{$bpchr}{$bp} = 0;
		}
		++$chrbpstatus{$bpchr}{$bp};
#		print "chrisstatus{$ischr}{$is}: $chrisstatus{$ischr}{$is}\n";
#		print "chrbpstatus{$bpchr}{$bp}: $chrbpstatus{$bpchr}{$bp}\n";
		my $isbp = "$ischr,$is,$bpchr,$bp";
		push @isbps, $isbp;
		$isbpline{$isbp} = $line;
		$isbppasscutoff{$isbp} = $passcutoff;
	}	
}
close CSV;

open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
print OUT "$header,Quality_control\n";
foreach my $isbp (@isbps) {
	my $quality = "N";
	my ($ischr, $is, $bpchr, $bp) = split /,/, $isbp;
	if ($chrisstatus{$ischr}{$is} > 1 or $chrbpstatus{$bpchr}{$bp} > 1) {
		$quality = "Y";
	}elsif ($isbppasscutoff{$isbp} < 3) {
		$quality = "Y";
	}
	print OUT "$isbpline{$isbp},$quality\n";
}
close OUT;
