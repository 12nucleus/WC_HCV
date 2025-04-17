#!/usr/bin/perl
use strict;


my $working_directory=$ARGV[0];
my $filename1=$ARGV[1];
my $output_dir = $ARGV[2];
my $bbmerge = "/usr/local/bin/bbmap/bbmerge.sh";
my $jellyfish = "/usr/local/bin/jellyfish";
my $cutadapt = "/usr/local/bin/cutadapt";


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
#system ("java -jar $working_directory/trimmomatic-0.36.jar PE -phred33 $filename1 $filename2 $prefix\_trimmed_R1_paired.fastq.gz R1_unpaired.fastq.gz $prefix\_trimmed_R2_paired.fastq.gz $working_directory/R2_unpaired.fastq.gz SLIDINGWINDOW:10:20\n");
system ("$cutadapt -G GGATATGATGATGAACTGGT -g ATGTGCCAGCTGCCGTTGGTGT -g GGATATGATGATGAACTGGT -G ATGTGCCAGCTGCCGTTGGTGT -o $output_dir/$prefix.1.fastq -p $output_dir/$prefix.2.fastq $filename1 $filename2 > $output_dir/cutadapt.log");
 


=for comment
This  block is to clean and merge the set of HCV reads as overlapping reads
with a minimum overlapp of 150bp and no missmatch allowed.
I do not trim the non-overlapped bases.  May lead to introduction of some sequencing errors
in the reads but should not cause any major issue down the line.  May need to be revisited 
if issue arises.
=cut	

### Combine reads with bbmerge, no missmatch allowed
system ("$bbmerge pfilter=1 minoverlap0=150 minoverlap=150 in1=$output_dir/$prefix.1.fastq in2=$output_dir/$prefix.2.fastq out=$working_directory/combined_reads/$prefix.fq.gz ziplevel=6 2> $output_dir/bbmerge.log");

print "Total reads\tPercent merged\t# HVR1\t# Rejected\tStatus\n";
open (FILE,"$output_dir/bbmerge.log");
while (my $in_log_file=<FILE>) {
	if ($in_log_file =~ 'Pairs') {
		chomp $in_log_file;
		my @split_log_line = split ('\t',$in_log_file);
		print "$split_log_line[1]\t";
		$in_log_file=<FILE>;
		chomp $in_log_file;
		@split_log_line = split ('\t',$in_log_file);
		print "$split_log_line[2]\t";	
	}
}
close FILE;

system ("rm $output_dir/$prefix.1.fastq");
system ("rm $output_dir/$prefix.2.fastq");
system ("rm $output_dir/bbmerge.log");
system ("rm $output_dir/cutadapt.log");
=for comment
This block is to consolidate sets of identical reads into a single sequence
as well as to verify that each consolidated reads is an HVR1 sequence (based on a dictionary
of HVR1 kmer set).  There is also a minimum cutoff ($n_minimum) for how many times a sequence is represented
in the dataset.
=cut


my %seq=();
my $total_reads=0;
#my $rejected_wrong_size=0;



my $read_file = "$working_directory/combined_reads/$prefix.fq.gz";
open (FILE,"gunzip -c $read_file |");


while (my $infile=<FILE>) {  ## read in gzip combined reads files
	$infile=<FILE>;
	++$total_reads;
	#if ((length($infile) >= 270) && (length($infile) <= 340)) { ## read lengths needs to be <270..340>  ## removed since subsequent filtering will catch any non-HVR1
		#my $truncated_infile = substr ($infile,11,999); ## remove the last 11nt and first 17nt corresponding to the primer locations
		#$truncated_infile = substr ($truncated_infile,0,-17);
		++$seq{$infile};	
	#} else {
	#	++$rejected_wrong_size;
	#}
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
 
my $n_minimum=15;
my @tmp_good_reads;
my $total=0;## total number of reads belonging to HVR1
my $rejected_minimum_representation=0;
my $rejected_kmer_matches=0;
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
			
		} else {
		$rejected_kmer_matches += $seq{$key};
		
		}
	} else {
		$rejected_minimum_representation += $seq{$key};
	}
++$tag;
}

close FILE;
print "$total\t$rejected_kmer_matches\t";

if ($total >= 50) { ## if total number of reads passing in the dataset > 500
	open (OUT,">$working_directory/HCV_fasta/$prefix\.fasta");
	print OUT @tmp_good_reads;	
	close OUT;	
	print "PASS\n";
	## create Kmer set from fasta file
	system ("$jellyfish count -m 125 -s 100M -t 1 $working_directory/HCV_fasta/$prefix\.fasta -o $output_dir/mer_counts.jf");
	system ("$jellyfish dump $output_dir/mer_counts.jf > $working_directory/HCV_Kmers/$prefix\.kmers");
	system ("gzip $working_directory/HCV_Kmers/$prefix\.kmers --force");
	system ("gzip $working_directory/HCV_fasta/$prefix\.fasta --force");
	system ("rm $output_dir/mer_counts.jf");	
} else {
	## report in statistics file a fail
	print "FAIL\n";
}





