#! /usr/bin/perl

use utf8;
use strict;
use warnings;
use File::Basename;
my $pwd = dirname(__FILE__);
die "perl $0 <*.csv> <exon.fa>\n" unless(@ARGV eq 2);
my %falength = get_fa_length($ARGV[1]);
my %genename;
my %cdsloc;
my %exons;
my %strands;
my %chrs;
open FA,$ARGV[1] or die $!;
while(<FA>){
	if(/^>/ && $_ =~ />(.*?)\s.*?CDS=(.*?)\s.*?loc:(.*?)\|.*?exons:(.*?)\s.*?segs:(.*)$/){
		my ($tran,$CDS,$chr,$exon,$seg) = $_ =~ />(.*?)\s.*?CDS=(.*?)\s.*?loc:(.*?)\|.*?exons:(.*?)\s.*?segs:(.*)$/;
		my $parent;
		if($_ =~ /gene=(.*?)\s/){
			$parent = $1;
		}else{
			$parent = $tran;
		}
		$genename{$tran} = $parent;
		$cdsloc{$tran} = $CDS;
		$chrs{$tran} = $chr;
		my $pos = index($_,"|",0);
		$pos = index($_,"|",$pos+1);
		my $strand = substr($_,$pos+1,1);
		$exons{$tran} = $exon;
		$strands{$tran} = $strand;
	}
}
sub get_raw_loc{
	my $loc = $_[0];
	my $strand = $_[1];
	my $seg = $_[2];
	my @segs = split(",",$seg);
	my @seg_length;
	my $summary_length = 0;
	my $raw_loc;
	if($strand eq "+"){
		for(my $i=0;$i<@segs;$i++){
			my($start,$end) = split("-",$segs[$i]);
			my $left_length = $loc - $summary_length;
			$summary_length = $summary_length + ($end-$start+1);
			if($summary_length >= $loc){
				$raw_loc = $start + $left_length-1;
				last;
			}
		}
	}elsif($strand eq "-"){
		for(my $i=scalar(@segs)-1;$i>=0;$i--){
			my($start,$end) = split("-",$segs[$i]);
			my $left_length = $loc - $summary_length;
			$summary_length = $summary_length + ($end-$start+1);
			if($summary_length >= $loc){
				$raw_loc = $end - $left_length + 1;
				last;
			}
		}
	}else{
		$raw_loc = "erro strand";
	}
	return $raw_loc;
}
open FILE,$ARGV[0] or die "failed to open file1";
<FILE>;
print "Transcript_ID\tGene_name\tTranscript_length\tcds_locs\tpause_codon\tcodon_pos\tchr\tchr_sit\treads_mapped\tpause_score\tcoverage\tz-score\tupstream_sequence\tdownstream_sequence\n";
while(<FILE>){
	chomp;
	my @as = split(",",$_);
	print $as[0]."\t".$genename{$as[0]}."\t".$falength{$as[0]}."\t".$cdsloc{$as[0]}."\t".$as[8]."\t".$as[1]."\t".$chrs{$as[0]}."\t".&get_raw_loc($as[1],$strands{$as[0]},$exons{$as[0]})."\t".$as[2]."\t".$as[3]."\t".$as[4]."\t".$as[7]."\t".$as[5]."\t".$as[6]."\n";
};
