#!/usr/bin/perl -w

##########################################################################################
# From mapping list by parsing _.sam file and R1_tid.txt file, mapping info for  
# each template
# Author: Wenjie Deng
# Date: 2017-12-12
########################################################################################## 

use strict;
use Getopt::Long;
use File::Basename;

my $usage = "usage: getTemplateMapping.pl hxb2MappingList R1tidReadnamefile R2tidReadnamefile outputTemplateMapping\n";

my $isfile = shift or die $usage;
my $tidfile = shift or die $usage;
my $r2tidfile = shift or die $usage;
my $outfile = shift or die $usage;

my (%namestatus, %nameIS, %nameRef, %nameDir, %r2tidcount);
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

open R2TID, $r2tidfile or die "couldn't open $r2tidfile: $!\n";
while (my $line = <R2TID>) {
	chomp $line;
	next if ($line =~ /^\s*$/ or $line =~ /^tid/);
	my @fields = split /\t/, $line;
	my $tid = $fields[0];
	my $counts = $fields[1];
	$r2tidcount{$tid} = $counts;
}
close R2TID;

open TID, $tidfile or die "couldn't open $tidfile: $!\n";
open OUT, ">", $outfile or die "couldn't open $outfile: $!\n";
print OUT "templateID,R2_count,R1_count,maping_count\n";
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
			}			
		}		
		print OUT "$tid,$r2tidcount{$tid},$counts,$actualcount\n";	
	}	
}
close TID;
close OUT;

