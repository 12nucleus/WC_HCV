@file=glob("/local/ngstb/NGSarchive-2015/150604_M01698_0164_000000000-AB83K-HCV-Flu-Dengue-wolfgang/HCV/BL*R1*.fastq.gz");

foreach $filename1(@file) {


	$filename2 = $filename1;
	$filename2 =~ s/\_R1\_/\_R2\_/g;  ## Second pair filename
	@split = split ('\/',$filename1);
	$total = @split;
	@split = split ('\_',$split[$total-1]);
	$prefix = $split[0];

	#system ("/local-homes/bioinformatics/erica/bbmap/bbmerge.sh pfilter=1 minoverlap0=225 minoverlap=225 in1=$filename1 in2=$filename2 out=/home/pxl10/Projects/PHGC/parker/combined_reads_old/$prefix.fq");

	system ("flash -t 5 -m 150 -M 500 -x 0 -o $prefix -d /home/pxl10/Projects/PHGC/parker/combined_reads $filename1 $filename2");
	system ("rm /home/pxl10/Projects/PHGC/parker/combined_reads/*.hist*");
	system ("rm /home/pxl10/Projects/PHGC/parker/combined_reads/*.notCombined*");

	}
	
