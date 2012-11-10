#!/usr/bin/perl -w

$path_to_db = "/home/ivanam/databases/ncbi_tax";

(-e $path_to_db) || die "$path_to_db not found\n";

chdir $path_to_db;

foreach $db  ("ncbi_taxonomy") {


    print "************************\n";
    print $db, "\n";
    chdir "$path_to_db/$db";

    print "making db  ... \n";
    $cmd = "mysqladmin -u root create $db";
    (system $cmd) && warn "error running $cmd\n";

    $sqlname = "$db.sql";
    (-e $sqlname) || die "$sqlname not found in $path_to_db/$db";
    $cmd = "mysql -u root $db < $sqlname";
    (system $cmd) && warn "error running $cmd\n";
=pod
mysqlimport loads tables from text files in various formats.  The base name of the
text file must be the name of the table that should be used.
=cut
    $cmd = "mysqlimport -u root ".
	"--fields_terminated_by=\'\\t|\\t\'  ".
	"--lines_terminated_by=\'\\t|\\n'\  ".
	"$db -L *.dmp";
    (system $cmd) && die "error running $cmd\n";

    print "\n";	
}
