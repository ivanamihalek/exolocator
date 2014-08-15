#! /usr/bin/perl -w 

$release_num = 74;

 $local_repository = 
    "/mnt/ensembl-mirror/release-$release_num/fasta";

chdir $local_repository;

@farm =  split "\n", `ls -d *_*`;

$ct = 0;
foreach $animal (@farm) {
    $ct ++;
    print "$ct   $animal \n";
}
print "----------------------------\n";
print " ... formating ... \n\n";

foreach $animal (@farm) {

    print $animal, "\n";


    foreach $dir  ( "pep",  "dna" ){

	chdir $local_repository;
	chdir "$animal/$dir";
	
	@fastas = split "\n",  `ls *.fa`;
	for $fasta (@fastas) {
	    print "\t $fasta\n";
	    $cmd = "formatdb -i $fasta -o T";
	    ( $dir eq "dna") && ( $cmd .= " -p F");
	    (system $cmd) &&
		die  " error running $cmd\n";
	}
    }
}
   
