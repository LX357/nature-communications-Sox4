#! /usr/bin/perl

use utf8;
use strict;
use warnings;
use Time::HiRes qw /usleep/;
die "perl $0 <sam> <*_rib_loc.txt> <output_prefix>\n" unless(@ARGV eq 3);
if($ARGV[0] =~ /.sam$/){
	open SAM,$ARGV[0] or die "failed to open $ARGV[0]";
}elsif($ARGV[0] =~ /.bam$/){
	open SAM,"samtools view -h $ARGV[0]|" or die $!;
}
open BAM,"|samtools view -Sb - > $ARGV[2].mark.bam" or die "failed to write $ARGV[2].mark.bam";
open BAMU,"|samtools view -Sb - > $ARGV[2].uniq.bam" or die "failed to write $ARGV[2].uniq.bam";
my $bedtools = `which bedtools`;
chomp($bedtools);
my $flag= "";
my $count = 0;
my $record = "";
my $first_line;
my $all=0;
while(<SAM>){
	if(/^@/){
		print BAM $_;
		print BAMU $_;
	}else{
		my @ars = split("\t",$_);
		if($ars[0] ne $flag){
			if($count <= 8 && $count != 0){
				#add NH:i:$count marks
				my @tmp = split("\n",$record);
				foreach my $value(@tmp){
					print BAM $value."\tNH:i:".$count."\n";
				}
				if($count == 1){
					chomp($first_line);
					print BAMU $first_line."\tNH:i:1\n";
					$all++;
				}
				$record = "";
			}
			$count = 1;
			$flag = $ars[0];
			$record = $_;
			$first_line = $_;
		}else{
			$record = $record.$_;
			$count++;
		}
	}
}
my @tmp = split("\n",$record);
foreach my $value(@tmp){
	print BAM $value."\tNH:i:".$count."\n";
}
if($count == 1){
	chomp($first_line);
	print BAMU $first_line."\tNH:i:1\n";
	$all++;
}
close BAMU;
close BAM;
`$bedtools intersect -a $ARGV[1] -b $ARGV[2].uniq.bam -F 1 -c > $ARGV[2].count.txt`;
open COUNT,"$ARGV[2].count.txt" or die  "failed to read $ARGV[2].count.txt";
my $coding = 0;
my $UTR3 = 0;
my $UTR5 = 0;
my $intron = 0;
while(<COUNT>){
	chomp;
	my @ars = split("\t",$_);
	if($ars[3] =~ /coding/){
		$coding += $ars[-1];
	}elsif($ars[3] =~ /5UTR/){
		$UTR5 += $ars[-1];
	}elsif($ars[3] =~ /3UTR/){
		$UTR3 += $ars[-1];
	}elsif($ars[3] =~ /intron/){
		$intron += $ars[-1];
	}
}
close COUNT;
my $other = $all - $coding - $UTR3 - $UTR5 - $intron;
my $pre = $ARGV[2];
$pre =~ s/.*\///g;
open TEMP, ">$ARGV[2].draw.txt" or die $!;
print TEMP "types\tnums
Coding	$coding
5UTR	$UTR5
3UTR	$UTR3
Intron	$intron
";
close TEMP;
