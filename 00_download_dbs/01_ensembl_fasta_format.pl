#! /usr/bin/perl -w 

 $local_repository = 
    "/mnt/ensembl-mirror/release-69/fasta";

chdir $local_repository;

@farm =  split "\n", `ls -d *_*`;

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
   
