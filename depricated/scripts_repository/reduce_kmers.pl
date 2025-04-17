use List::Util 'shuffle';

@files = glob ("*.kmers.gz");
$nfiles=@files;
$x=0;
foreach $file(@files) {
	print "$file\n";
	$prefix1=$file;
	$prefix1 =~ s/\.kmers\.gz//g;
	open(FILE, "gunzip -c $file |");
	@infile = <FILE>;
	close FILE;
	$join = join  ('',@infile);
	@split = split ('\>',$join);
	 @shuffled_indexes = shuffle(0..$#split);
	 @pick_indexes = @shuffled_indexes[ 0 .. 99999 ];
	 open (OUT,">$prefix1.reduced.kmers");
	 print OUT join '>',@split[ @pick_indexes ];
	 close OUT;
	 
}