#!/usr/bin/perl -w
# ftp://ftp.ensembl.org/pub/current_mysql/ncbi_taxonomy/
# note that table dumps (ncbi_taxa_name.txt.gz, ncbi_taxa_node.txt.gz)
# are separate from the schema (ncbi_taxonomy.sql.gz)
use strict;
use warnings;

my $version = 101;

my $path_to_db = "/storage/databases/ensembl-$version/ncbi_taxonomy";
(-e $path_to_db) || die "$path_to_db not found\n";

chdir $path_to_db;



foreach my $db  ("ncbi_taxonomy") {


    print "************************\n";
    print $db, "\n";
    chdir "$path_to_db";

    print "making db  ... \n";
    # do not check for error here (if the database exists)"
    my $cmd  = "mysqladmin  --login-path=tcga  drop $db --force";
    $cmd  = "mysqladmin  --login-path=tcga create $db";
 
    print "\t $cmd\n";
    (system $cmd) && warn "error running $cmd - moving on\n\n";


    my $sqlname = "$db\_$version.sql";
    (-e $sqlname) || die "$sqlname not found in $path_to_db";
    # as of mysql v 8,  NFORMATION_SCHEMA.SESSION_VARIABLES became performance_schema.session_variables
    `sed 's/INFORMATION_SCHEMA.SESSION_VARIABLES/performance_schema.session_variables/g' -i  $sqlname `;
    $cmd = "mysql  --login-path=tcga   $db < $sqlname";
    (system $cmd) && die "error running $cmd - exiting\n";

    #mysqlimport loads tables from text files in various formats.  The base name of the
    #text file must be the name of the table that should be used.

    $cmd = "mysqlimport --login-path=tcga  ".
	"$db -L *.txt";
    (system $cmd) && die "error running $cmd\n";

    print "\n";	
}
