#!/usr/bin/perl

use strict;

(my $indir = shift @ARGV) or die $!;

print join("\t", "PF", "total", "duplicates", "mapped", "prop_paired", "secondary",
	"singletons", "chrM", "Qlt20", "FLlt38", "FLgt2K",
	"mean_fsize", "nuc_free", "mono_nuc", "di_nuc", "tri_nuc",
	"peaks", "frip", "min_plen", "median_plen", "mean_plen", "max_plen"), "\n";

my %fs = &read_flagstat("${indir}/Aligned.out.mdup.flagstat");
my %is = &read_idxstats("${indir}/Aligned.out.mdup.idxstats");
my %fl = &read_fsize("${indir}/Aligned.out.filt.fs");
my %p = &read_peaks("${indir}/broad_peaks.count", "${indir}/broad_peaks.count.summary");

my $fassigned = $p{unassigned} ?  $p{assigned}/($p{assigned}+$p{unassigned}) : 0;
print join("\t", $fl{reads},
	$fs{total}, $fs{duplicates}, $fs{mapped}, $fs{'properly paired'},
	$fs{secondary}, $fs{singletons}, $is{chrM},
	$fs{total} - &cat("${indir}/Aligned.out.mdup.q20"),
	&cat("${indir}/Aligned.out.mdup.lt38"),
	&cat("${indir}/Aligned.out.mdup.gt2k"),
	$fl{mean_len}, $fl{nucleosome_free}, $fl{mono_nucleosome},
	$fl{di_nucleosome}, $fl{tri_nucleosome},
	$p{peaks}, $fassigned,
	$p{min}, $p{median}, $p{mean}, $p{max}), "\n";

sub read_flagstat {
	my $file = shift;
	my %res = ();
	open IN, "< $file" or die $!;
	while(<IN>) {
		s/\s+$//;
		s/ in / /;
		s/\s\(.+$//;
		/^(\d+) \+ 0 (.+)$/ || die;
		$res{$2} = $1;
	}
	close IN;
	return %res;
}

sub read_idxstats {
	my $file = shift;
	my %res = ();
	open IN, "< $file" or die $!;
	while(<IN>) {
		my($name, $n) = (split /\t/)[0,2];
		$res{$name} = $n;
	}
	close IN;
	return %res;
}

sub read_fsize {
	my $file = shift;
	open IN, "< $file" or die $!;
	$_ = <IN>; s/\s+$//;
	my @cols = split /\t/;
	$_ = <IN>; s/\s+$//;
	my @vals = split /\t/;
	close IN;
	my %res = ();
	foreach my $c (@cols) {
		my $v = shift @vals;
		$res{$c} = $v;
	}
	return %res;
}

sub read_peaks {
	my($count, $summary) = @_;
	my %res = (peaks=>0, min=>0, median=>0, mean=>0, max=>0,
	           assigned=>0, unassigned=>0);
	return %res unless -f $count;
	my $sum = 0;
	my @lens = ();
	open IN, "< $count" or die $!;
	<IN>; <IN>;
	$res{peaks} = 0;
	while(<IN>) {
		my $l = (split /\t/)[5];
		$sum += $l;
		push @lens, $l;
	}
	close IN;
	my @sorted = sort {$a <=> $b} @lens;
	$res{peaks} = scalar(@sorted);
	$res{min} = $sorted[0];
	$res{median} = $sorted[int(scalar(@sorted)/2)];
	$res{mean} = $sum / scalar(@sorted);
	$res{max} = $sorted[$#sorted];

	open IN, "< $summary" or die $!;
	while(<IN>) {
		s/\s+$//;
		if(/^Assigned\s+(\d+)$/) {
			$res{assigned} = $1;
		} elsif(/^Unassigned_NoFeatures\s(\d+)$/) {
			$res{unassigned} = $1;
		}
	}
	close IN;
	return %res;
}

sub cat {
	my $file = shift;
	open IN, "< $file" or die $!." ".$file;
	my $res = <IN>; $res =~ s/\s+$//;
	close IN;
	return $res;
}
