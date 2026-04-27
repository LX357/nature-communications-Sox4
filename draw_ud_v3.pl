#! /usr/bin/perl

use utf8;
use strict;
use warnings;
use lib "/home/guojianhui/lib";
use Fzuguo;
use File::Basename;
my $pwd = dirname(__FILE__);
die "perl $0 <dir/*.bam> <gtf>\n" unless(@ARGV eq 2);
my $basename = basename($ARGV[0]);
my $path = $ARGV[0];
$path =~ s/$basename$//;
$basename =~ s/\.bam//g;
`mkdir -p $path/$basename`;
`ln -s $ARGV[0] $path/$basename/$basename.bam`;
print ("Rscript do.r $basename $path $ARGV[1]");
`Rscript do.r $basename $path $ARGV[1]`;

