#!/usr/bin/perl

use strict;
use IPC::Run 'run'; ##run [ "command", "arguments", "here" ], ">", \my $stdout;
use File::Temp qw/tempdir/;

=for comment
usage :

HCV_Master_script_V0.1.pl Path_to_reads

=cut




my $working_directory="/Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline";
my @file = glob ("$ARGV[0]*");
my $email = $ARGV[1];
$email =~ s/\_\_at\_\_/\@/g;
$email =~ s/ /\,/g;
my $prefix;
my @samples;
my $combine_reads_script = "$working_directory/HCV_combine_reads_V0.1.pl";
my $kmer_distancer_script = "$working_directory/HCV_kmers_distancer_V0.1.pl";
my $test_transmission = "$working_directory/HCV_transmission_test_V0.1.pl";
my $rundir;

=for comments
this block of code is to prepare the reads, combine overlapping reads
assess HCV content and produce a kmer sketch
=cut
my @files_to_send;
foreach my $filename1(@file) {
	my $filename2;
	if (($filename1 =~ '\_R1\_001\.fastq\.gz') && ($filename1 !~ 'Undeter|undeter')) {
		$rundir = tempdir ( 'run_XXXXXXXXXXXX',DIR => $working_directory, TMPDIR=> 1,CLEANUP => 0);  ## set and create a random running directory 
		$filename2 = $filename1;
		$filename2 =~ s/\_R1\_/\_R2\_/g;		
		my @split_path = split ('\/',$filename1);
		my @split_filename = split ('\_',$split_path[-1]);
		$prefix = $split_filename[0];
		open (OUT,">$working_directory/Reports/$prefix\_report.out");
		push (@files_to_send,"$working_directory/Reports/$prefix\_report.out");
		print OUT "Processing sample $prefix:\n\n";
		push (@samples,$prefix);				
		## Call Read Combiner/corrector script ##
		
		run [ "$combine_reads_script", "$working_directory", "$filename1" ,"$rundir"], ">", \my $stdout;
		print OUT "$stdout\n\n";
		close OUT;
		system ("rm -r $rundir");
	}
}

=for comments
this block of code is compare different kerm sketch and determine if there
are any overlapping of kmer set to determine mixture between samples
and assume transmission took place
=cut

### run kmer distancer only when all new dataset are processed

foreach my $to_process(@samples) {
	open (OUT,">>$working_directory/Reports/$to_process\_report.out");
	
	if (-e "$working_directory/HCV_Kmers/$to_process.kmers.gz") {
		run [ "$kmer_distancer_script", "$working_directory", "$to_process" ], ">", \my $stdout;
		print OUT "$stdout\n\n";
		my @split_kmer_stats = split ('\n',$stdout);
		my @check_transmission=(); 
		my $potential_links=0;
		foreach my $line_kmer_stats(@split_kmer_stats) {
			my @split_line_kmer = split ('\t', $line_kmer_stats);
			if ($split_line_kmer[5] >= 0.05) {
				++$potential_links;
				push (@check_transmission,$split_line_kmer[1]);
			}
		}
		if ($potential_links == 0) {
			print OUT "\nNo Potential transmission cases detected for sample $to_process\n";
		} else {
			print OUT "\n$potential_links potential transmission case(s) detected for sample $to_process\n\n";
			foreach my $input1(@check_transmission) {
				print OUT "Testing $to_process VS $input1\n--------------------------------\n\n";
				run [ "$test_transmission", "/Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline$to_process.fasta.gz", "/Volumes/DOH_HOME/pxl10/Projects/HCV_pipeline$input1.fasta.gz", "$working_directory" ],">", \my $stdout;
				print OUT "$stdout\n\n";
				push (@files_to_send,"$working_directory/Reports/$prefix.pdf");
			}
		}
		
		
		
		
	} else { 
		print OUT "cannot find $working_directory/HCV_Kmers/$to_process.kmers.gz\n"; 
	}
	close OUT;
}

my $join_file_list = join " ",@files_to_send;
exit;
#system ("cd $working_directory/Reports/; zip Reports.zip -j $join_file_list");
#system ("sendemail -t '$email' -m \"Please do not reply to this email.  For any questions or issues, contact the pipeline administrator at pascal.lapierre\@health.ny.gov\" -u \"VirCapSeq pipeline reports $date\" -a $working_directory/Reports/Reports.zip -f bioinfo\@health.ny.gov");
#system ("rm $working_directory/Reports/Reports.zip");

