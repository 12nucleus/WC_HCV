#!/usr/bin/perl -w
use List::Util 'shuffle';


@files = glob ("*.fasta");
foreach $file(@files) {
	$reduce=15000000;
	$prefix = $file;
	$prefix =~ s/\.fasta//g;
	system ("jellyfish count -m 125 -s 100M -t 10 $file");
	system ("jellyfish dump mer_counts.jf > $prefix.kmers");
	system ("rm mer_counts.jf");
	
	#open(FILE, "$prefix.kmers");
	#@infile = <FILE>;
	#close FILE;
	#$join = join  ('',@infile);
	#@split = split ('\>',$join);
	#$n = @split;
	#if ($n < $reduce) {
	#	open (OUT,">$prefix.reduced.kmers");
	#	print OUT $join;
	#} else {
		
		
		
	#@shuffled_indexes = shuffle(0..$#split);
	#@pick_indexes = @shuffled_indexes[ 0 .. $reduce-1 ];
	#open (OUT,">$prefix.reduced.kmers");
	#print OUT ">";
	#print OUT join '>',@split[ @pick_indexes ];
	#}
	#close OUT;
	
	
	
	
	
	
	
	
	#system ("gzip $prefix.kmers -f");
}



@files = glob ("*.kmers");
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