#!/usr/bin/perl -w
use List::Util 'shuffle';


@files = glob ("/local/ngs/illumina_results/170505_M01698_0013_000000000-B4BKY-HCV-Vimi/Data/Intensities/BaseCalls/GP*_R1_*.fastq.gz");
foreach $file(@files) {
	@split = split ('\/',$file);
	@split = split ('\.',$split[-1]);
	$prefix = $split[0];
	$file2 = $file;
	$file2 =~ s/\_R1\_/\_R2\_/g;
	system ("zcat $file $file2 > tmp.fq");
	
	system ("jellyfish count -m 125 -s 100M -t 10 tmp.fq");
	system ("jellyfish dump mer_counts.jf > $prefix.kkmers");
	system ("rm mer_counts.jf");
	system ("rm tmp.fq");

	
	
	exit;
	
	
	#system ("gzip $prefix.kmers -f");
}



@files = glob ("*.kkmers");
$nfiles=@files;
$x=0;
do {
	%set1=();$total1=0;
	$prefix1=$files[$x];
	$prefix1 =~ s/\.kmers//g;
	
	#open(FILE, "gunzip -c $files[$x] |");
	open (FILE,$files[$x]);
	while ($infile=<FILE>) {
	
		if ($infile !~ '>') {
			chomp $infile;
			$set1{$infile}=1;
			++$total1;
		}
	}	
	close FILE;	
	
	$y = $x+1;
	do {
		
		$total2=0;
		%set2=();
		$prefix2=$files[$y];
		$prefix2 =~ s/\.kmers//g;
	
		
		#open(FILE, "gunzip -c $files[$y] |");
		open (FILE,$files[$y]);
		while ($infile=<FILE>) {
			if ($infile !~ '>') {
				chomp $infile;
				$set2{$infile}=1;
				++$total2;
			}
		}	
		close FILE;
		
		

		$overlapp=0;
		foreach ( keys %set1 ) {
			if ( exists $set2{$_} ) {
				++$overlapp;
			}
	
		}
		if ($overlapp >= 1) {
			
			print "$prefix1\t$prefix2\t$overlapp\t$total1\t$total2\n";
		}
		++$y;
	}
	until $y==$nfiles;
++$x;
}
until $x==$nfiles-1;