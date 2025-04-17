#!/usr/bin/perl -w
use List::Util 'shuffle';
#use strict;

##compare prefix to all other kmer files
## create a figure of scores in R?
## if 1 or more other kmer above thresholds, retain or run inter-intra script;

my $working_directory=$ARGV[0];
my $query=$ARGV[1];
my @kmer_files = glob ("$working_directory/HCV_Kmers/*.kmers.gz");

my %set1=();
my $total1=0;
my $prefix1=$query;
print "S1\tS2\tKmer Overlaps\tTotal S1\tTotal S2\tRatio\n";

open(FILE, "gunzip -c $working_directory/HCV_Kmers/$query.kmers.gz |");
while (my $infile=<FILE>) {
	if ($infile !~ '>') {
		chomp $infile;
		$set1{$infile}=1;
		++$total1;
	}
}	
close FILE;




foreach my $subject(@kmer_files) {
	my $total2=0;
	my %set2=();
	my $prefix2;
	my @split_directory = split ('\/',$subject);
	my @split_filename = split ('\.',$split_directory[-1]);
	$prefix2 = $split_filename[0];
		if ($prefix1 !~ $prefix2) {
		open(FILE, "gunzip -c $subject |");
	
		while (my $infile=<FILE>) {
			if ($infile !~ '>') {
				chomp $infile;
				$set2{$infile}=1;
				++$total2;
			}
		}	
		close FILE;

		my $overlapp=0;
		foreach ( keys %set1 ) {
			#print "$set1{$_}\n";
			if ( exists $set2{$_} ) {			
				++$overlapp;
			}
		}
	
		if ($overlapp >= 1) {
			my $r;
			if ($total1 < $total2) {
				$r = sprintf ("%.5f",$overlapp / $total1);
			} else {
				$r = sprintf ("%.5f",$overlapp / $total2);
			}

			print "$prefix1\t$prefix2\t$overlapp\t$total1\t$total2\t$r";
			if ($r >= 0.05) {
				print "*";
			}
			print "\n";
			
		}
	}
}
print "* indicated Kmer overlapp above threshold of 0.05\n";