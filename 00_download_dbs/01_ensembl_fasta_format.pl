#! /usr/bin/perl -w 
# not sure if this is needed any more at all - at some point I switched to

use warnings;
use strict;
use Cwd;
my $release_num = 110;

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

	#$animal eq "danio_rerio" && next;

	print "\n----------------------------\n";
	print " ... formating  $animal ... \n\n";

    foreach my $dir  ( "dna", 'pep' ){

		chdir $local_repository;
		chdir "$animal/$dir";
        print ("I am in ", getcwd, "\n");
		my @fastas = split "\n",  `ls *.fa.gz`;
		for my $fasta (@fastas) {
		    my $new_db_name = $fasta;
		    $new_db_name =~ s/.gz$//;
		    my $title = $new_db_name;
			my $cmd;
			if ( $dir eq "dna") {
				if ((-e "$new_db_name.nsq" && ! -z  "$new_db_name.nsq") || (-e "$new_db_name.00.nsq" && ! -z  "$new_db_name.00.nsq")) {
					print "\t found formatted $fasta \n";

				} else {
					print "\t formatting $fasta\n";
					# $cmd = "$db_formatting_tool -in $fasta -dbtype nucl -parse_seqids";
					$cmd = "gunzip -c $fasta | $db_formatting_tool -out $new_db_name -dbtype nucl -parse_seqids";
					$cmd .= " -title $title";
					(system $cmd) && die " error running $cmd\n";
				}

			} else {
				if (-e "$new_db_name.psq" && ! -z  "$new_db_name.psq") {
					print "\t found formatted $fasta \n";

				} else {
					print "\t formatting $fasta\n";
					#$cmd = "$db_formatting_tool -in $fasta -dbtype prot -parse_seqids";
					$cmd = "gunzip -c $fasta | $db_formatting_tool -out $new_db_name  -dbtype prot -parse_seqids";
					$cmd .= " -title $title";
					(system $cmd) && die " error running $cmd\n";
				}
			}
		}
    }
}
