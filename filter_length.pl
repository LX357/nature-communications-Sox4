#! /usr/bin/perl

use utf8;
use strict;
use warnings;
use File::Basename qw/basename dirname/;
use File::Spec::Functions qw/rel2abs/;
die "perl $0 <fq.gz> <minlenth> <maxlength> <output_prefix>\n" unless(@ARGV eq 4);
my $path = dirname(rel2abs($ARGV[0]));
open FILE,"gzip -dc $ARGV[0]|" or die "failed to open $ARGV[0]";
open OUT,"|gzip > $path/$ARGV[3]_ft.fq.gz" or die  "failed to write $ARGV[3]";
my %count;
my $allcount;
while(<FILE>){
	my $l1 = $_;
	my $l2 = <FILE>;
	my $l3 = <FILE>;
	my $l4 = <FILE>;
	my $length = length($l2)-1;
	if($length>=$ARGV[1] && $length<=$ARGV[2]){
		$count{$length}++;
		$allcount++;
		print OUT $l1.$l2.$l3.$l4;
	}
}
open TEMP,"> $path/$ARGV[3].RF_length.txt";
print TEMP "length\tpercent\n";
foreach my $key(sort {$a<=>$b} keys %count){
	my $percent = sprintf("%.2f",100*$count{$key}/$allcount);
	print TEMP $key."\t$percent\n";
}
close TEMP;
