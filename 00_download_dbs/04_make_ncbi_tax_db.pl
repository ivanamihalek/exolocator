#!/usr/bin/perl -w
# ftp://ftp.ensembl.org/pub/current_mysql/ncbi_taxonomy/
use strict;
use warnings;
my $path_to_db = "/storage/databases/ensembl-97/ncbi_tax";
(-e $path_to_db) || die "$path_to_db not found\n";

chdir $path_to_db;



foreach my $db  ("ncbi_taxonomy") {


    print "************************\n";
    print $db, "\n";
    chdir "$path_to_db/$db";

    print "making db  ... \n";
    my $cmd  = "mysqladmin  --login-path=client create $db";
 
    print "\t $cmd\n";
    (system $cmd) && warn "error running $cmd\n";


    my $sqlname = "$db.sql";
    (-e $sqlname) || die "$sqlname not found in $path_to_db/$db";
    $cmd = "mysql  --login-path=client   $db < $sqlname";
    (system $cmd) && warn "error running $cmd\n";

    #mysqlimport loads tables from text files in various formats.  The base name of the
    #text file must be the name of the table that should be used.

    $cmd = "mysqlimport --login-path=ivana  ".
	"--fields_terminated_by=\'\\t|\\t\'  ".
	"--lines_terminated_by=\'\\t|\\n'\  ".
	"$db -L *.dmp";
    (system $cmd) && die "error running $cmd\n";

    print "\n";	
}
