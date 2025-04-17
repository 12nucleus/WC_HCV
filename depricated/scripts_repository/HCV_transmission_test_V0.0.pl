#!/usr/bin/perl;
#use strict;
use Statistics::TTest;
#open (OUT,">test_results.out");
use Statistics::Descriptive;
my $stat = Statistics::Descriptive::Full->new();
use Math::CDF;
open (OUT,">clost.out");
@clear;


my $file1=$ARGV[0];
my $file2=$ARGV[1];
#$file1="/home/pxl10/Projects/PHGC/parker/kmer_approach/old_testfiles/BL-HCVNGS014.fasta";
#$file2="/home/pxl10/Projects/PHGC/parker/kmer_approach/old_testfiles/GP00002.fasta";


@tmp_split = split ('\/',$file1);
@tmp_split = split ('\.',$tmp_split[-1]);
$sample1 = $tmp_split[0];
@tmp_split = split ('\/',$file2);
@tmp_split = split ('\.',$tmp_split[-1]);
$sample2 = $tmp_split[0];


system ("zcat $file1 > tmp1.fa");
system ("zcat $file2 > tmp2.fa");



system ("/usr/local/bin/mafft --adjustdirection --auto --quiet --thread 1 --reorder tmp1.fa > data1.out 2> log.out");
system ("/usr/local/bin/mafft --adjustdirection --auto --quiet --thread 1 --reorder tmp2.fa > data2.out 2> log.out");
system ("cat tmp1.fa tmp2.fa > data3.fasta");
system ("/usr/local/bin/mafft --adjustdirection --auto --quiet --thread 1 --reorder data3.fasta > data3.out 2> log.out");


$data1 = "data1.out";
$data2 = "data2.out";
$data3 = "data3.out"; 

my $ham_sample = $sample1;
hamming($data1);
$stat = Statistics::Descriptive::Full->new();
$stat->add_data(@clear);
$stat->add_data(@array);
@array_data1=@array; 

$mean_set1 = $stat->mean();
$stdev_set1 = $stat->standard_deviation();
$set1_count = @array;

#print "$sample1:\t$mean_set1\t$stdev_set1\n";
@head_set1=@header;
#my %set1=('count' =>$set1_count,'mean' =>$mean_set1,'variance' =>$stdev_set1);
#print "$set1_count,$mean_set1,$stdev_set1\n";
my $ham_sample = $sample2;
hamming($data2);
$stat = Statistics::Descriptive::Full->new();
$stat->add_data(@clear);
$stat->add_data(@array); 
@array_data2=@array; 


#print join "\n",@array;
#exit;
$mean_set2 = $stat->mean();
$stdev_set2 = $stat->standard_deviation();
#print "$sample2:\t$mean_set2\t$stdev_set2\n\n";
$set2_count = @array;
@head_set2=@header;
#exit;
hamming($data3);

#= for comments
#p-values for Z-scores:
#z of 1.65 -> 0.1
#z of 1.96 -> 0.05
#z of 2.58 -> 0.01
#z of 0.8416 -> 0.2 (80 percentile cutoff ingroup, i.e, z < 0.8416 only if you belong into the 80 percentile,
# you are flagged as belonging to this group..hash cutoff)
#=cut
$z_cut=1.282;
$z_cut=0.8416;
## compare set 1 versus set 2
$oo=0;
$oi=0;
$io=0;
$ii=0;
foreach $seq1(@head_set1) {
	@array1=();
	@avg=();
	#$seq1 = "408_GP00004_24";
	$q1=0;
	$q2=0;
	
	foreach $seq2(@head_set1) {
		
		if ($seq1 !~ $seq2) {
	
 				$z = (($hash{$seq1}{$seq2} - $mean_set1)/$stdev_set1);
 				$zz=$hash{$seq1}{$seq2};
				push (@avg,$zz);
 		
 		#if ($z >= 2.58) {
 		#	$q1=1;
 		#print "$seq1 $seq2\n";
 		#exit;
 		$prob = Math::CDF::pnorm($z);
 		 #if ($z < 0) {
 		#	$prob = 1-$prob;
 		#}
 		#print "Q1: $seq1 $seq2\t$z\t$prob\t$hash{$seq1}{$seq2}\n";
 		#}

 		#exit;
 		#print "$seq1\t$seq2\t$z\t$hash{$seq1}{$seq2}\n";
 		#print "$seq1\t$seq2\t$hash{$seq1}{$seq2}\n";
 		push (@array1,$z);
 		
 		#	push (@array1,$hash{$seq1}{$seq2});	
 		}
 	}
 	#$array1_count=@array1;
 	
 	
 	$stat = Statistics::Descriptive::Full->new();
 	$stat->add_data(@clear);
	$stat->add_data(@array1); 
	$mean_sample1 = $stat->mean();
	if ($mean_sample1 >= $z_cut) {
		$q1=1;
	}
	$stat = Statistics::Descriptive::Full->new();
 	$stat->add_data(@clear);
	$stat->add_data(@avg);
	$mean_sample1 = $stat->mean();
	
	

	#$stdev_sample1 = $stat->standard_deviation();
 	#my %sample1=('count' =>$array1_count,'mean' =>$mean_sample1, 'variance' =>$stdev_sample1);
 	#print "$array1_count\t$mean_sample1\t$stdev_sample1\n";
	
 	#my $ttest = new Statistics::TTest;  
	#$ttest->load_data(\@array_data1,\@array1);  
	#print "$seq1 ", $ttest->{t_prob},"\n\n\n";
	
 	
 	@array2=();
 	@avg=();
 	foreach $seq2(@head_set2) {
		
		if ($seq1 !~ $seq2) {

 			$z = (($hash{$seq1}{$seq2} - $mean_set2)/$stdev_set2);
 			 				$zz=$hash{$seq1}{$seq2};
				push (@avg,$zz);
 		#if ($z >= 2.58) {
 		#	$q2=1;
 		$prob = Math::CDF::pnorm($z);

 		#print "Q2: $seq1 $seq2\t$z\t$prob\t$hash{$seq1}{$seq2}\n";
 		#}
 		#print "$seq1\t$seq2\t$z\t$hash{$seq1}{$seq2}\n";
 		push (@array2,$z);	
 		#push (@array2,$hash{$seq1}{$seq2});
 			#print "$seq1\t$seq2\t$hash{$seq1}{$seq2}\n";
 		}
 	}
 	
 	$stat = Statistics::Descriptive::Full->new();
 	$stat->add_data(@clear);
	$stat->add_data(@array2); 
	$mean_sample2 = $stat->mean();
	#$stdev_sample2 = $stat->standard_deviation();
	#$array2_count=@array2;
	
	if ($mean_sample2 >= $z_cut) {
		$q2=1;
	}
	
	#if ($mean_sample1 < $mean_sample2) {
	#	$q1=0;$q2=1;
	#} else {
	#	$q1=1;$q2=0;
	#}
		$stat = Statistics::Descriptive::Full->new();
 	$stat->add_data(@clear);
	$stat->add_data(@avg);
	$mean_sample1 = $stat->mean();
	
 	#my %sample2=('count' =>$array2_count,'mean' =>$mean_sample2, 'variance' =>$stdev_sample2);
 	#print "$set2_count,$mean_set2,$stdev_set2\n";
 	#print "$array2_count\t$mean_sample2\t$stdev_sample2\n";
 	#print join "\n",@array2;
 	
 	#$ttest2 = new Statistics::TTest; 
	#$ttest2->load_data(\@array_data2,\@array2);  
	#print $ttest2->{t_prob},"\n";
 	#print "$seq1\t$q1\t$q2\n";
 	if (($q1==0) && ($q2==0)) {
 		++$oo;
 	} elsif (($q1==1) && ($q2==0)) {
 		++$io;
 	} elsif (($q1==0) && ($q2==1)) {
 		++$oi;	
 	} else {
 		++$ii;
 	}
		
}	

print "$sample1:\t$oo\t$oi\t$io\t$ii\n";

$tot = $oo+$oi+$io+$ii;
$result_s1="Failed";
if ($oo/$tot >= 0.8) {
	$result_s1="Mixed";
} elsif ($oi/$tot >= 0.5) {
	$result_s1="Homogenious";
} elsif ($io/$tot >= 0.5) {
	$result_s1="recipient";
} elsif (($ii > $oo) && ($ii > $oi) && ($ii > $io)) {
	$result_s1="too many unclassified reads";
} 













 	
 ## compare set 2 versus set 1
$oo=0;
$oi=0;
$io=0;
$ii=0;
foreach $seq1(@head_set2) {
	@array1=();
	#$seq1 = "637_GP00004_17";
	$q1=0;
	$q2=0;
	foreach $seq2(@head_set1) {
		
		if ($seq1 !~ $seq2) {
 		$z = (($hash{$seq1}{$seq2} - $mean_set1)/$stdev_set1);
 		
 		#if ($z >= 2.58) {
 		#	$q1=1;
 		
 		$prob = Math::CDF::pnorm($z);

 		#print "Q1: $seq1 $seq2\t$z\t$prob\t$hash{$seq1}{$seq2}\n";
 		#}

 		#exit;
 		#print "$seq1\t$seq2\t$z\t$hash{$seq1}{$seq2}\n";
 		#print "$seq1\t$seq2\t$hash{$seq1}{$seq2}\n";
 		push (@array1,$z);	
 		#push (@array1,$hash{$seq1}{$seq2});	
 		}
 	}
 	#$array1_count=@array1;
 	
 	
 	$stat = Statistics::Descriptive::Full->new();
 	$stat->add_data(@clear);
	$stat->add_data(@array1); 
	$mean_sample1 = $stat->mean();
	if ($mean_sample1 >= $z_cut) {
		$q1=1;
	}

	
	
	#$stdev_sample1 = $stat->standard_deviation();
 	#my %sample1=('count' =>$array1_count,'mean' =>$mean_sample1, 'variance' =>$stdev_sample1);
 	#print "$array1_count\t$mean_sample1\t$stdev_sample1\n";
	
 	#my $ttest = new Statistics::TTest;  
	#$ttest->load_data(\@array_data1,\@array1);  
	#print "$seq1 ", $ttest->{t_prob},"\n\n\n";

 	
 	@array2=();
 	foreach $seq2(@head_set2) {
		
		if ($seq1 !~ $seq2) {
 		$z = (($hash{$seq1}{$seq2} - $mean_set2)/$stdev_set2);
 		#if ($z >= 2.58) {
 		#	$q2=1;
 		$prob = Math::CDF::pnorm($z);

 		
 		#print "Q2: $seq1 $seq2\t$z\t$prob\t$hash{$seq1}{$seq2}\n";
 		#}
 		#print "$seq1\t$seq2\t$z\t$hash{$seq1}{$seq2}\n";
 		push (@array2,$z);	
 		#push (@array2,$hash{$seq1}{$seq2});
 			#print "$seq1\t$seq2\t$hash{$seq1}{$seq2}\n";
 		}
 	}
 	
 	$stat = Statistics::Descriptive::Full->new();
 	$stat->add_data(@clear);
	$stat->add_data(@array2); 
	$mean_sample2 = $stat->mean();
	#$stdev_sample2 = $stat->standard_deviation();
	#$array2_count=@array2;
	
	if ($mean_sample2 >= $z_cut) {
		$q2=1;
	}

	
	#if ($mean_sample1 < $mean_sample2) {
#		$q1=0;$q2=1;
#	} else {
#		$q1=1;$q2=0;
#	}
 	#my %sample2=('count' =>$array2_count,'mean' =>$mean_sample2, 'variance' =>$stdev_sample2);
 	#print "$set2_count,$mean_set2,$stdev_set2\n";
 	#print "$array2_count\t$mean_sample2\t$stdev_sample2\n";
 	#print join "\n",@array2;
 	
 	#$ttest2 = new Statistics::TTest; 
	#$ttest2->load_data(\@array_data2,\@array2);  
	#print $ttest2->{t_prob},"\n";
 	#print "$seq1\t$q1\t$q2\n";
 	if (($q1==0) && ($q2==0)) {
 		++$oo;
 	} elsif (($q1==1) && ($q2==0)) {
 		++$io;
 	} elsif (($q1==0) && ($q2==1)) {
 		++$oi;	
 	} else {
 		++$ii;
 	}	
 	
}


print "$sample2:\t$oo\t$oi\t$io\t$ii\n";

$tot = $oo+$oi+$io+$ii;
$result_s2="Failed";
if ($oo/$tot >= 0.8) {
	$result_s2="Mixed";
} elsif ($oi/$tot >= 0.5) {
	$result_s2="recipient";
} elsif ($io/$tot >= 0.5) {
	$result_s2="Homogenious";
} elsif (($ii > $oo) && ($ii > $oi) && ($ii > $io)) {
	$result_s2="too many unclassified reads";
} 

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
print "Who am I? God?  I cant decide..\n";
print OUT "Who am I? God?  I cant decide..\n";

}
print "______________________________________\n\n";
close OUT;
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
			#exit;
			
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



