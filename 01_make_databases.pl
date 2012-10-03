#!/usr/bin/perl -w


$path_to_db = "/home/ivanam/databases/ensembl/release-68";


chdir $path_to_db;


@dbs = split "\n", `ls -d [n-z]*`;



foreach $db (@dbs) {

    print "************************\n";
    print $db, "\n";
    chdir "$path_to_db/$db";
    print "unzipppig ... \n";
    `gunzip *.gz`;

    print "making db  ... \n";
    $cmd = "mysqladmin -u root create $db";
    (system $cmd) && warn "error running $cmd\n";

    $cmd = "mysql -u root $db < $db.sql";
    (system $cmd) && warn "error running $cmd\n";

    $cmd = "mysqlimport -u root --fields_escaped_by=\\\\ $db -L *.txt";
    (system $cmd) && die "error running $cmd\n";




    print "\n";	
}
