#! /usr/bin/perl -w 

$release_num = 91;

#$local_repository = "/mnt/ensembl-mirror/release-$release_num/fasta";
$local_repository    = "/data/ensembl-$release_num/fasta";
$db_formatting_tool  = "/usr/bin/makeblastdb";
# makeblastdb allows for Maximum file size: 1000000000B
# note this is not enough for 'toplevel' files
for ( $local_repository,  $db_formatting_tool) {
    ( -e $_) ||  die "$_ not found\n";
}

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
	    if ( $dir eq "dna") {
		    $cmd = "$db_formatting_tool -in $fasta -dbtype nucl -parse_seqids";
	    } else {
		    $cmd = "$db_formatting_tool -in $fasta -dbtype prot -parse_seqids";
	    }
	    (system $cmd) &&
		die  " error running $cmd\n";
	}
    }
}
