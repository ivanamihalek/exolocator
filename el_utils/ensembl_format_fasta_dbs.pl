#! /usr/bin/perl -w 

 $local_repository = 
    "/home/ivanam/databases/ensembl/fasta";

chdir $local_repository;

@farm =  split "\n", `ls -d *_*`;

foreach $animal (@farm) {

    print $animal, "\n";


    #foreach $dir  ( "pep",  "dna" ){
    foreach  $dir  ( "pep" ){

	chdir $local_repository;
	chdir "$animal/$dir";
    
	@fastas = split "\n",  `ls *.fa`;
	for $fasta (@fastas) {
	    print "\t $fasta\n";
	    $cmd = "formatdb -i $fasta -o T\n";
	    (system $cmd) &&
		die  " error running $cmd\n";
	}
    }
}
   
