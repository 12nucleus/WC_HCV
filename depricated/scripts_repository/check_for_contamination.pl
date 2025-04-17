open (FILE,"kmer_ref_counts.out");
while ($infile=<FILE>) {
if ($infile !~ '>') {
	chomp $infile;
	push (@kmer_check,$infile);
	}
}
close FILE;


@glob = glob ("GP00003.fasta");
foreach $file(@glob) {

	open (FILE,$file);
	print "$file\t";
	$low=1;
	$high=0;
	$avg=0;
	$for_avg=0;
	$n=0;
	while ($infile=<FILE>) {
		if ($infile !~ '>') {
			++$n;
			$match=0;
			$total=0;
			foreach $check(@kmer_check) {
				if ($infile =~ $check) {
					++$match;
					}
				++$total;
				}
			$r = $match/$total;
			if ($r < 0.6) {
				print ">$n\_$r\n$infile\n";
				}
			$for_avg += $r;
			
			if ($low > $r) {
				$low=$r;
			}
			if ($high < $r) {
				$high=$r;
			}
			
			
			
		}
	}
	close FILE;
	exit;
	$avg = $for_avg/$n;
	print "$low\t$high\t$avg\n";
	
	exit;
	
	
}