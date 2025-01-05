#!/usr/bin/perl

use strict;

@ARGV || die $0." BAMs";
my @bams = ();
foreach (@ARGV) {
	push @bams, $_;
}

my $sum = 0;
my $n = 0;
my $n_free = 0;
my $mono_n = 0;
my $di_n = 0;
my $tri_n = 0;

foreach my $bam (@bams) {
	open IN, "samtools view $bam | cut -f9 |" or die $!;
	while(my $f = <IN>) {
		$f =~ s/\s+$//;
		$f = abs $f;
		$sum += $f;
		$n ++;
		if($f >= 100 && $f <= 180) {
			$n_free ++;
		} elsif($f >= 180 && $f <= 247) {
			$mono_n ++;
		} elsif($f >= 315 && $f <= 473) {
			$di_n ++;
		} elsif($f >= 558 && $f <= 615) {
			$tri_n ++;
		}
	}
	close IN;
}

print join("\t", "reads", "mean_len", "nucleosome_free", "mono_nucleosome", "di_nucleosome",
	"tri_nucleosome"), "\n";
print join("\t", $n, $sum / $n, $n_free, $mono_n, $di_n, $tri_n), "\n";

exit;
