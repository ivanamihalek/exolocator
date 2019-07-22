#! /usr/bin/perl -w 

use warnings;
use strict;
my $release_num = 97;

my $local_repository    = "/storage/databases/ensembl-$release_num/fasta";
my $db_formatting_tool  = "/usr/bin/makeblastdb";
# makeblastdb allows for Maximum file size: 1000000000B
# note this is not enough for 'toplevel' files
for ( $local_repository,  $db_formatting_tool) {
    ( -e $_) ||  die "$_ not found\n";
}

chdir $local_repository;

my @farm =  split "\n", `ls -d *_*`;



foreach my $animal (@farm) {

	print "----------------------------\n";
	print " ... formating  $animal... \n\n";

    foreach my $dir  ( "pep",  "dna" ){

		chdir $local_repository;
		chdir "$animal/$dir";

		my @fastas = split "\n",  `ls *.fa`;
		for my $fasta (@fastas) {
			print "\t $fasta\n";
			my $cmd;
			if ( $dir eq "dna") {
				$cmd = "$db_formatting_tool -in $fasta -dbtype nucl -parse_seqids";
			} else {
				$cmd = "$db_formatting_tool -in $fasta -dbtype prot -parse_seqids";
			}
			(system $cmd) && die  " error running $cmd\n";
		}
    }
}
