use Digest::MurmurHash3 qw( murmur32 );
use Digest::xxHash64 qw(xxhash32 xxhash32_hex xxhash64 xxhash64_hex);
use Array::Utils qw(:all);



$seed=934595;


open (FILE,"GP00013.kmers");

while ($infile=<FILE>) {
	if ($infile !~ '>') {
		chomp $infile;
		
		$h=xxhash32( $infile,$seed );
		#print "$h\n";
		push (@array1,$h);
		
	}
}
close FILE;
#@sorted1 = sort { $a <=> $b } @array1;
#$t1=@sorted2;


open (FILE,"GP00011.kmers");
while ($infile=<FILE>) {
	if ($infile !~ '>') {
		chomp $infile;
		
		$h=xxhash32( $infile,$seed );
		#print "$h\n";
		push (@array2,$h);
		
	}
}
close FILE;

$nn = intersect(@array1, @array2);
#$nn=@isect;
print $nn;
open (OUT,">test.out");
print OUT join "\n",@array1;

exit;

@sorted2 = sort { $a <=> $b } @array2;
$t2=@sorted2;


if ($t1 < $t2) {
	$samples=$t1;
} else {
	$samples=$t2;
}

if ($samples < $total) {
	$total = $samples;
}

$x=0;
do {
	$hash1{$sorted1[$x]} = 1;
	#print "$sorted[$x]\n";
	
	++$x;
}
until $x==$total;




$x=0;
do {
	$hash2{$sorted2[$x]} = 1;
	#print "$sorted[$x]\n";
	
	++$x;
}
until $x==$total;

foreach my $key (keys(%hash1)) {
	if (exists $hash2{$key}) {
		++$shared;
	}
}

print "$shared\n";