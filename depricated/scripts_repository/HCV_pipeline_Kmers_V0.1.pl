#!/usr/bin/perl
use strict;



my $filename1=$ARGV[0];
my $working_directory="/home/pxl10/Projects/PHGC/parker/kmer_approach";


=for comment
This small loop loads in @kmer_check a general set of 15 mers derived from a small
portion of the HVR1 region coming from a general sets of HCV strains. this set of kmer 
will serve to identified if we are dealing with an actual HVR1 sequence or non-specific 
contaminant.
=cut


my @kmer_check;
open (FILE,"$working_directory/kmer_ref_counts.out");
while (my $infile=<FILE>) {
	if ($infile !~ '>') {
		chomp $infile;
		push (@kmer_check,$infile);
	}
}
close FILE;


my $filename2 = $filename1;
$filename2 =~ s/\_R1\_/\_R2\_/g;  ## Second pair filename
my @split_directory = split ('\/',$filename1);
my @split_filename = split ('\_',$split_directory[-1]);
my $prefix = $split_filename[0];

#system ("/local-homes/bioinformatics/erica/bbmap/bbmerge.sh pfilter=1 minoverlap0=225 minoverlap=225 in1=$filename1 in2=$filename2 out=/home/pxl10/Projects/PHGC/parker/combined_reads_old/$prefix.fq");


=for comment
This  block is to clean and merge the set of HCV reads as overlapping reads
with a minimum overlapp of 150bp and no missmatch allowed.
I do not trim the non-overlapped bases.  May lead to introduction of some sequencing errors
in the reads but should not cause any major issue down the line.  May need to be revisited 
if issue arises.
=cut	

### Combine reads with flash
system ("flash -t 5 -m 150 -M 500 -x 0 -z -o $prefix -d $working_directory/combined_reads $filename1 $filename2");
system ("rm $working_directory/combined_reads/*.hist*");
system ("rm $working_directory/combined_reads/*.notCombined*");


=for comment
This block is to consolidate sets of identical reads into a single sequence
as well as to verify that each consolidated reads is an HVR1 sequence (based on a dictionary
of HVR1 kmer set).  There is also a minimum cutoff ($n_minimum) for how many times a sequence is represented
in the dataset.
=cut

my $n_minimum=15;
my %seq=();

my $rejected=0;

open (OUT,">$working_directory/HCV_fasta/$prefix\.fasta");	
my $read_file = "$working_directory/combined_reads/$prefix.extendedFrags.fastq.gz";
open (FILE,"gunzip -c $read_file |");


while (my $infile=<FILE>) {  ## read in gzip combined reads files
	$infile=<FILE>;

	if ((length($infile) >= 270) && (length($infile) <= 340)) { ## read lengths needs to be <270..340>
		my $truncated_infile = substr ($infile,11,999); ## remove the last 11nt and first 17nt corresponding to the primer locations
		$truncated_infile = substr ($truncated_infile,0,-17);
		++$seq{$truncated_infile};		
	} else {
		++$rejected;
	}
	$infile=<FILE>;
	$infile=<FILE>;
} 
close FILE;
my $tag=1; ## unique tag to be added to the read name in order to avoid any name collision


=for comment
this last loop sort the hash of sequence reads from the biggest count to the lowest
and check through a series of kmers derived from a diverse population of HCH HVR1 sequences
to determined if the reads are actual HVR1 sequences instead of contamination or non-specific
amplicons.  To pass, the read needs to have at least $n_minimum representatives and 15 matches to
the HVR1 kmer database
=cut
 

my @tmp_good_reads;
my $total=0;## total number of reads belonging to HVR1
foreach my $key (sort { $seq{$b} <=> $seq{$a} } keys %seq) { ## not necessary to sort but just look tidier in the file
	if (($seq{$key} >= $n_minimum)) {
		my $match=0;  ## number of kmer match to our kmer reference database		  
		foreach my $check(@kmer_check) {
			my $rcheck = $check;
			$rcheck = reverse $rcheck;
			$rcheck =~ tr/ATCG/TAGC/;
			if (($key =~ $check) || ($key =~ $rcheck)) { ## check for both the forward and reverse kmer
				++$match;
			}
		}		
		if ($match >= 15) {  ## need at least 15 kmer match 
			$total += $seq{$key};
			push (@tmp_good_reads,">$tag\_$prefix\_$seq{$key}\n$key\n");
		}
	}
++$tag;
}




if ($total >= 500) { ## if total number of reads passing in the dataset > 500
	print OUT @tmp_good_reads;
} else {
	## report in statistics file a fail
}

close FILE;
close OUT;

