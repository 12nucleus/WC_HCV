
open (FILE,$ARGV[0]);
$in=<FILE>;
$string=$in;
while ($in=<FILE>) {
if ($in=~ '>') {
$string =~ s/\-/X/g;
$string .= "\n";
push (@con,$string);
$string=$in;
} else {
chomp $in;
$string .=$in;
}
}
$string .= "\n";
$string =~ s/\-/X/g;

push (@con,$string);
#@in=<FILE>;
close FILE;

$num = @con;


open (OUT,">dist.matrix");
open (OUTT,">ident.pmatrix");

open (OUTTT,">dist_normalized.matrix");
open (OUTTTT,">ident_normalized.pmatrix");


foreach $e(@con) {
@s1 = split ('\n',$e);
print OUT "\t$s1[0]";
print OUTT "\t$s1[0]";
print OUTTT "\t$s1[0]";
print OUTTTT "\t$s1[0]";

}

print OUTTT "\n";
print OUTTTT "\n";
print OUT "\n";
print OUTT "\n";
$x=0;

	local $now = time;
	### get current date ###

do {


$y=0;
@s1 = split ('\n',$con[$x]);
print "$s1[0]\t";
print OUT "$s1[0]";
print OUTT "$s1[0]";
print OUTTT "$s1[0]";
print OUTTTT "$s1[0]";
$ln = length($s1[1]);

do {

if ($y==$x) {
print OUT "\t0";
print OUTT "\t0";
print OUTTT "\t0";
print OUTTTT "\t0";
} elsif ($y < $x) {
print OUT "\t$dd[$y][$x]";
print OUTT "\t$ppp[$y][$x]";
print OUTTT "\t$dd_norm[$y][$x]";
print OUTTTT "\t$ppp[$y][$x]";

} else {
@s2 = split ('\n',$con[$y]);
print "$s2[0]\n";
$same=hd($s1[1],$s2[1]);
$nn=hdn($s1[1],$s2[1]);
$n1=n($s1[1]);
$n2=n($s2[1]);
$end = $same - (($n1-$nn) + ($n2-$nn));
$v = valid($s1[1],$s2[1]);
$pp = sprintf ("%.6f",$end/($end + $v));
$total_length = length($s1[1]);
$normal_factor = ($end+$v)/$total_length;
$norm_end = sprintf ("%.6f", $end/$normal_factor);




print OUT "\t$end";
print OUTT "\t$pp";

print OUTTT "\t$norm_end";
print OUTTTT "\t$pp";
$dd_norm[$x][$y] = $norm_end;
$dd[$x][$y] = $end;
$ppp[$x][$y] = $pp;


}
++$y;
}
until $y>=$num;

print OUT "\n";
print OUTT "\n";
print OUTTT "\n";
print OUTTTT "\n";







++$x;
}
until $x>=$num;
$now = time - $now;
printf ("\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));



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
    return ($_[0] & $_[1]) =~ tr/ATCGX//;
}