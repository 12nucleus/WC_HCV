#!/usr/bin/perl
use strict;
use File::Temp qw/tempdir/;

my $file1=$ARGV[0];
my $file2=$ARGV[1];
my $working_directory=$ARGV[2];
my $rundir = tempdir ( 'run_XXXXXXXXXXXX',DIR => $working_directory, TMPDIR=> 1,CLEANUP => 0);  ## set and create a random running directory

open (OUT,">$rundir/test_results.out");
open (LEGEND,">$rundir/legend.txt");
use Statistics::Descriptive;
use Math::CDF;


print LEGEND "Samples\tCat\tSize\n";

### extract sample names
my @tmp_split_filename = split ('\/',$file1);
@tmp_split_filename = split ('\.',$tmp_split_filename[-1]);
my $sample1 = $tmp_split_filename[0];
@tmp_split_filename = split ('\/',$file2);
@tmp_split_filename = split ('\.',$tmp_split_filename[-1]);
my $sample2 = $tmp_split_filename[0];
my (%hash,$mean_sample1,$mean_sample2,$seq2,$z,$tot,$result_s1,$result_s2,@con,@header,$in,$string,$num,$x,$y,@s1,@s2,$n1,$n2,$nn,$v,$pp,$same,$end,$total_length,$ln,@ppp);
## gather read fasta files for dataset1 and dataset2 and individually align
system ("zcat $file1 > $rundir/tmp1.fa");
system ("zcat $file2 > $rundir/tmp2.fa");
system ("/usr/local/bin/mafft --adjustdirection --auto --quiet --thread 1 --reorder $rundir/tmp1.fa > $rundir/data1.out 2> $rundir/mafft_log.out");
system ("/usr/local/bin/mafft --adjustdirection --auto --quiet --thread 1 --reorder $rundir/tmp2.fa > $rundir/data2.out 2>> $rundir/mafft_log.out");

## merge dataset1 and dataset2 and align
system ("cat $rundir/tmp1.fa $rundir/tmp2.fa > $rundir/data3.fasta");
system ("/usr/local/bin/mafft --adjustdirection --auto --quiet --thread 1 --reorder $rundir/data3.fasta > $rundir/data3.out 2>> $rundir/log.out");

my $dataset1 = "$rundir/data1.out";
my $dataset2 = "$rundir/data2.out";
my $dataset_all = "$rundir/data3.out"; 

=for comments
these next two sets of commands is use to calculate all the pairwise distance (percent identity) within each dataset.
These distances are kept in an array (@array_data) as well as the name of each read analyzed (@head_dataset).
A mean and standard deviation that will be used to get z-probabilities is calculated from that same array.
The hamming subroutine is storing a hash using the pair of header name as the key, and the distance as the value for later retrieval
=cut
my @array;
hamming($dataset1); 
my $stat = Statistics::Descriptive::Full->new();
$stat->add_data(@array);
my @array_data1=@array; 
my $mean_dataset1 = $stat->mean();
my $stdev_dataset1 = $stat->standard_deviation();
my $dataset1_count = @array;
my @head_dataset1=@header;


hamming($dataset2);
$stat = Statistics::Descriptive::Full->new();
$stat->add_data(@array); 
my @array_dataset2=@array; 
my $mean_dataset2 = $stat->mean();
my $stdev_dataset2 = $stat->standard_deviation();
my $dataset2_count = @array;
my @head_dataset2=@header;


## subroutine to get hash of all pairwise comparisons
hamming($dataset_all);

=for comments
	p-values for Z-scores:
	z of 1.65 -> 0.1
	z of 1.96 -> 0.05
	z of 2.58 -> 0.01
	z of 0.8416 -> 0.2 (80 percentile cutoff ingroup, i.e, z < 0.8416 only if you belong into the 80 percentile,
 	you are flagged as belonging to this group..hash cutoff)
=cut

my $z_cut=0.8416;  ## cutoff to determine if a read is part of a sample or not.  This is very harsh cutoff since if a read
## falls in the upper 20 percentile, we say we cannot be sure if he is part of the group (20% chance of not being part of this group)


## compare set 1 versus set 2
my $AA=0;  #AA = read Belong to group 1 and 2, i.e., intermixture of samples
my $AB=0;  #AB = read Belong to group 1 but not to group 2
my $BA=0;  #BA = read Belong to group 2 but not to group 1
my $BB=0;  #BB = read do not cluster in any of the two groups, i.e., outlier
	my $t_dist=0;
	my $avg_dist=0;
	my $total_data=0;
foreach my $seq1(@head_dataset1) { # iterate through header for dataset1
	
	my @name_split = split ('\_',$seq1);
	my $scaled_size = log($name_split[2])/log(10);
	print LEGEND "$seq1\t2\t$scaled_size\n";  ## generate legend file for the PCA plot with seq name,2(color) and bubble size
	
	
	$t_dist=0;
	$avg_dist=0;
	$total_data=0;

	my @zscore_dataset1=(); # array to store all z-score for the read analyzed
	my $q1=0;
	my $q2=0;
	foreach my $seq2(@head_dataset1) {	# iterate through header for dataset1 for pairwise comparison within dataset1
		if ($seq1 !~ $seq2) {
 			my $z = (($hash{$seq1}{$seq2} - $mean_dataset1)/$stdev_dataset1);  ## calculate the Zscore
 			#my $prob = Math::CDF::pnorm($z);
 			push (@zscore_dataset1,$z);
 			++$total_data;
 			$t_dist+= $hash{$seq1}{$seq2};
 		}
 	}


	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@zscore_dataset1); 
	$mean_sample1 = $stat->mean();
	if ($mean_sample1 >= $z_cut) {
		$q1=1;
	}
	
	$avg_dist = $t_dist/$total_data;
	#print OUT_DIST "$seq1\t1\t$avg_dist\n";
	$t_dist=0;
	$avg_dist=0;
	$total_data=0;	
	
	my @zscore_dataset2=();	
 	foreach $seq2(@head_dataset2) {
		if ($seq1 !~ $seq2) {
 			$z = (($hash{$seq1}{$seq2} - $mean_dataset2)/$stdev_dataset2);
 			#$prob = Math::CDF::pnorm($z);
 			push (@zscore_dataset2,$z);
 			 ++$total_data;
 			$t_dist+= $hash{$seq1}{$seq2};
 		}
 	}
 	$avg_dist = $t_dist/$total_data;
 	#print OUT_DIST "$seq1\t2\t$avg_dist\n";
	$t_dist=0;
	$avg_dist=0;
	$total_data=0;

	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@zscore_dataset2); 
	$mean_sample2 = $stat->mean();
	if ($mean_sample2 >= $z_cut) {
		$q2=1;
	}
	
 	if (($q1==0) && ($q2==0)) {
 		++$AA;
 	} elsif (($q1==1) && ($q2==0)) {
 		++$BA;
 	} elsif (($q1==0) && ($q2==1)) {
 		++$AB;	
 	} else {
 		++$BB;
 	}
		
}	

print "Sample\tIntermixture\tGroup1\tGroup2\tUnclustered\n--------------------------------------------\n";
print "$sample1:\t$AA\t$AB\t$BA\t$BB\n";

$tot = $AA+$AB+$BA+$BB;
$result_s1="Failed";
if ($AA/$tot >= 0.8) {
	$result_s1="Mixed";
} elsif ($AB/$tot >= 0.5) {
	$result_s1="Homogenious";
} elsif ($BA/$tot >= 0.5) {
	$result_s1="recipient";
} elsif (($BB > $AA) && ($BB > $AB) && ($BB > $BA)) {
	$result_s1="tAA many unclassified reads";
} 





## compare set 2 versus set 1
my $AA=0;  #AA = read Belong to group 1 and 2, i.e., intermixture of samples
my $AB=0;  #AB = read Belong to group 1 but not to group 2
my $BA=0;  #BA = read Belong to group 2 but not to group 1
my $BB=0;  #BB = read do not cluster in any of the two groups, i.e., outlier

foreach my $seq1(@head_dataset2) { # iterate through header for dataset2
	my @name_split = split ('\_',$seq1);
	my $scaled_size = log($name_split[2])/log(5);
	print LEGEND "$seq1\t4\t$scaled_size\n";
	
	$t_dist=0;
	$avg_dist=0;
	$total_data=0;
	my @zscore_dataset1=(); # array to store all z-score for the read analyzed
	my $q1=0;
	my $q2=0;
	foreach my $seq2(@head_dataset1) {	# iterate through header for dataset2 for pairwise comparison within dataset1
		if ($seq1 !~ $seq2) {
 			my $z = (($hash{$seq1}{$seq2} - $mean_dataset1)/$stdev_dataset1);
 			#my $prob = Math::CDF::pnorm($z);
 			push (@zscore_dataset1,$z);
 			 ++$total_data;
 			$t_dist+= $hash{$seq1}{$seq2};
 		}
 	}


	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@zscore_dataset1); 
	$mean_sample1 = $stat->mean();
	if ($mean_sample1 >= $z_cut) {
		$q2=1;
	}
	$avg_dist = $t_dist/$total_data;
	#print OUT_DIST "$seq1\t3\t$avg_dist\n";
	$t_dist=0;
	$avg_dist=0;
	$total_data=0;

	my @zscore_dataset2=();	
 	foreach $seq2(@head_dataset2) {
		if ($seq1 !~ $seq2) {
 			$z = (($hash{$seq1}{$seq2} - $mean_dataset2)/$stdev_dataset2);
 			#$prob = Math::CDF::pnorm($z);
 			push (@zscore_dataset2,$z);
 			 ++$total_data;
 			$t_dist+= $hash{$seq1}{$seq2};
 		}
 	}
	
	$stat = Statistics::Descriptive::Full->new();
	$stat->add_data(@zscore_dataset2); 
	$mean_sample2 = $stat->mean();
	if ($mean_sample2 >= $z_cut) {
		$q1=1;
	}
	$avg_dist = $t_dist/$total_data;
	#print OUT_DIST "$seq1\t4\t$avg_dist\n";
	$t_dist=0;
	$avg_dist=0;
	$total_data=0;
 	if (($q1==0) && ($q2==0)) {
 		++$AA;
 	} elsif (($q1==1) && ($q2==0)) {
 		++$AB;
 	} elsif (($q1==0) && ($q2==1)) {
 		++$BA;	
 	} else {
 		++$BB;
 	}
		
}	

print "$sample2:\t$AA\t$AB\t$BA\t$BB\n";

$tot = $AA+$AB+$BA+$BB;
$result_s2="Failed";  # set 'failed' as default state
if ($AA/$tot >= 0.2) {
	$result_s2="Mixed";  ##  if more that 20% of the total reads in mixed state
} elsif ($AB/$tot >= 0.5) {
	$result_s2="recipient"; ## if more that 50% matching 1-1
} elsif ($BA/$tot >= 0.5) {
	$result_s2="Homogenious";## if more that 50% matching 1-2
} elsif (($BB > $AA) && ($BB > $AB) && ($BB > $BA)) {
	$result_s2="tAA many unclassified reads"; # if too many reads in unmatched
} 

#close OUT_DIST;

close LEGEND;



print "\n\n$result_s1\t$result_s2\n\n";
print "$sample1 vs $sample2:\n";
print OUT "$sample1 vs $sample2:\t";
if (($result_s1 =~ 'Failed') || ($result_s2 =~ 'Failed')) {
	print "\nFailed to confirm data transmission\n";
	print OUT "Failed to confirm data transmission\n"
} elsif (($result_s1 =~ 'Mixed') && ($result_s2 =~ 'Mixed')) {
	print "\nComplete intersample mixture\n";
	print OUT "Complete intersample mixture\n";
} elsif (($result_s1 =~ 'Mixed') && ($result_s2 =~ 'Homogenious')) {
	print "\nTransmission from $sample2 to $sample1\n";
	print OUT "Transmission from $sample2 to $sample1\n";
} elsif (($result_s2 =~ 'Mixed') && ($result_s1 =~ 'Homogenious')) {
	print OUT "Transmission from $sample1 to $sample2\n";
	print "\nTransmission from $sample1 to $sample2\n";
} elsif (($result_s2 =~ 'Homogenious') && ($result_s1 =~ 'Homogenious')) {
	print OUT "No transmission detected\n";
	print "\nNo transmission detected\n";
} else {
print "Could Not Determine\n";
print OUT "Could Not Determine\n";

}
print "______________________________________\n\n";
close OUT;

close LEGEND;
system ("cd $rundir/; /usr/local/bin/hamming $rundir/data3.out");
open (R_OUT,">$rundir/rscript.txt");
print R_OUT "library(Rtsne)\n";
print R_OUT "library(ggplot2)\n";
print R_OUT "m <- as.matrix(read.table(file=\"$rundir/ident.pmatrix\", sep='\\t', header=T)[-1])\n";
#print R_OUT "fit <- cmdscale(m,eig=T,k=2)\n";
print R_OUT "fit <- Rtsne(m,dims=2, perplexity=30,check_duplicates = FALSE)\n";
#print R_OUT "y <- fit\$points\[,1\]\n";
#print R_OUT "y <- fit\$points\[,2\]\n";

print R_OUT "x <- fit\$Y\[,1\]\n";
print R_OUT "y <- fit\$Y\[,2\]\n";
print R_OUT "Table <- (read.table(file=\"$rundir/legend.txt\", sep='\\t', header=T))\n";
print R_OUT "df1 = data.frame(x=x, y=y)\n";
print R_OUT "pdf(file=\"$working_directory/Reports/Rtsne\_$sample1\_$sample2.pdf\",onefile=FALSE)\n";
print R_OUT "plot.new()\n";
print R_OUT "ggplot(df1, aes(x=x, y=y)) + geom_point(size=Table\$Size,color=\"black\",fill=Table\$Cat,shape=21,alpha=0.5)\n";
print R_OUT "legend(\"topright\",c(\"$sample1\",\"$sample2\"),pt\.bg=c(\"blue\",\"red\"),pch=c(21,21),pt.lwd=0.5,cex=0.5,bg=\"white\")\n";
print R_OUT "dev.off()\n";
close R_OUT;
system ("R CMD BATCH $rundir/rscript.txt");
system ("rm -r $rundir");




exit;


sub hamming {
	@array=(); 
	@con=();
	@header=();
	open (FILE,$_[0]);
	$in=<FILE>;
	$string=$in;
	$in =~ s/\>//g;
	$in =~ s/\_R\_//g;
	chomp $in;
	push (@header,$in);
	while ($in=<FILE>) {
		if ($in=~ '>') {
			
			$string .= "\n";
			push (@con,$string);
			$string=$in;
			chomp $in;
			$in =~ s/\>//g;
			$in =~ s/\_R\_//g;
			push (@header,$in);
		} else {
			chomp $in;
			$in = uc($in);
			$string .=$in;
			#$string =~ s/\-/X/g;
		}
	}
	$string .= "\n";
	#$string =~ s/\-/X/g;
	push (@con,$string);
	close FILE;
	$num = @con;


	$x=0;
	do {
		$y=$x+1;
		@s1 = split ('\n',$con[$x]);
		$s1[0] =~ s/\>//g;
		$s1[0] =~ s/\_R\_//g;
		$ln = length($s1[1]);

		do {

			$pp="";
			$end="";
			@s2 = split ('\n',$con[$y]);
			$s2[0] =~ s/\>//g;
			$s2[0] =~ s/\_R\_//g;
			$same=hd($s1[1],$s2[1]);
			$nn=hdn($s1[1],$s2[1]);
			$n1=n($s1[1]);
			$n2=n($s2[1]);
			$end = $same - (($n1-$nn) + ($n2-$nn));
			$v = valid($s1[1],$s2[1]);
			$pp = sprintf ("%.6f",$end/($end + $v));
			$total_length = length($s1[1]);

			push (@array,$pp);
			
			#print OUT "$ham_sample\t$pp\n";
			$hash{$s1[0]}{$s2[0]} = $pp;
			$hash{$s2[0]}{$s1[0]} = $pp;
			#print "$s1[0] $s2[0] $pp\n";

			$ppp[$x][$y] = $pp;
			++$y;
		}
		until $y>=$num;
   
		++$x;
	}
	until $x==$num-1;
}


sub hd {
    return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}
sub hdn {
    return ($_[0] & $_[1]) =~ tr/N//;
}
sub n {
    return ($_[0]) =~ tr/N//;
}
sub valid {
    return ($_[0] & $_[1]) =~ tr/ATCG\-//;
}


