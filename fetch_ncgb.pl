#!/usr/bin/perl
# fetch seq from ncgb

use strict;
use Getopt::Std;
use vars qw($opt_t $opt_s $opt_h);
getopts('t:s:h');

## program and options
my $types=$opt_t ? $opt_t : "";  # ncRNA types, concatenated by comma
my $org=$opt_s ? $opt_s : "";
my $help=$opt_h;

## USAGE 
my $usage = <<USAGE;
Program: Extract Seq From ncgb
Usage: perl $0 <options> ncgb.fa     
	-t [STR] types, eg: 'rRNA,tRNA,snRNA,snoRNA,scRNA',
			 default: ALL
	-s [STR] organism name, eg. 'Homo sapiens',
			 default: ALL
	-h       Help
Author: liqb <liqb\@genomics.org.cn> 2007-9-21
        miaoxin <xmiao\@genedenovo.com> 2015-11-06
USAGE

my $database=shift;
unless($database) {
	print $usage;
	exit;
}

$types=~s/\s//g;
$types=~tr/,/|/;
$org=~s/^\s*//; $org=~s/\s*$//;
$types = "" unless ($types);
$org = "" unless ($org);

#print "type: $types\n";
#print "org: $org\n";

open DB, $database || die $!;
my $flag=0;
while (my $line = <DB>) {
	if ($line =~ /^>(\S+)/) {
		my @wd = split(/_/, $1);
		my @context = split(/,/, $line);
		if ($wd[1]=~/$types/i && $context[0]=~/\b$org\b/) {
			$flag=1;
			print $line;
		} else {
			$flag = 0;
		}
	} elsif ($flag) {
		print $line;
	}
}

