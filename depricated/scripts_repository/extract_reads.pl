open (FILE,"kmer_ref_counts.out");
while ($infile=<FILE>) {
if ($infile !~ '>') {
	chomp $infile;
	push (@kmer_check,$infile);
	}
}
close FILE;


@file = glob ("/home/pxl10/Projects/PHGC/parker/combined_reads_old/GP00018*.extendedFrags.fastq"); 

$x=0;
open (LOG,">log.out");





foreach $ff(@file) {
@tmp = split ('\/',$ff);
@tmp = split ('\.',$tmp[-1]);
$seed = int(rand 99999);
#system ("seqtk sample -s seed $ff 50000 > subsampled.fq");

#system ("/home/pxl10/Projects/Ecoli_stools/bbmaps/bbduk.sh in=$ff outm=filtered.fq ref=HCV_ref.fa overwrite=true");

#system ("/usr/local/bin/bbmap/bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 path=/home/pxl10/Projects/PHGC/parker/kmer_approach/human_ref/ qtrim=rl trimq=10 untrim -Xmx23g in=$ff outu=clean.fq");# outm=human.fq


open (OUT,">$tmp[0]\.fasta");

$rejected=0;
%seq=();
@name = split ('\.',$ff);


@p = split ('\/',$name[0]);
$name[0] = $p[-1];

open (FILE,$ff);
$total=0;
while ($in=<FILE>) {
$in=<FILE>;
if ((length($in) >= 270) && (length($in) <= 340)) {
$inn = substr ($in,11,999);
$inn = substr ($inn,0,-17);
#$inn=$in;
++$seq{$inn};
++$total;
} else {
++$rejected;
}
$in=<FILE>;
$in=<FILE>;

} 

$x=1;


foreach my $key (sort { $seq{$b} <=> $seq{$a} } keys %seq) {
#print "$seq{$key}\n";
	#print "$seq{$key}\n";
	if (($seq{$key} >= 10)) {
#if (($x <= 10) && ($seq{$key} > 20)) {
	

		$match=0;
		$total=0;
		foreach $check(@kmer_check) {
			$rcheck = $check;
			$rcheck = reverse $rcheck;
			$rcheck =~ tr/ATCG/TAGC/;
			if (($key =~ $check) || ($key =~ $rcheck)) {
				++$match;
				#print "$x\_$name[0]\_$seq{$key}\t$check\t$rcheck\n";
			}
			++$total;
		}
		#$r = $match/$total;
#		print "$x\_$name[0]\_$seq{$key}\t$match\t$total\t$r\n";
		if ($match >= 15) {
		$m=$seq{$key};
		#$key .= $key;
		
		print OUT ">$x\_$name[0]\_$m\n$key\n";
		}
	}
++$x;
}



#foreach $key (keys %seq) {
#$pc = $seq{$key}/$total;
#if ($pc >= 0.0000) {
#if ($seq{$key} > 150) {

#print OUT ">$x\_$name[0]\_$seq{$key}\n$key\n";
#print "$seq{$key}\n";
#++$x;
#}
#}
close FILE;
#print LOG "$ff\tRejected = $rejected\n";
close OUT;

}


#system ("muscle -in data.fasta -out data.aln");


#system ("cat data.fasta HCV_hyper_controls.fasta > full_dataset.fasta");



