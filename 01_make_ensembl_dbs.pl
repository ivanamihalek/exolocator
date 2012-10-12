#!/usr/bin/perl -w
# one extra database that is needed,
# never made it a part of this script:
# ensembl_compara_xx, the tables needed are homology, homology_member, member and genome_db

$path_to_db = "/home/ivanam/databases/ensembl/mysql";

(-e $path_to_db) || die "$path_to_db not found\n";

chdir $path_to_db;


@dbs = split "\n", `ls -d *`;



foreach $db (@dbs) {

    print "************************\n";
    print $db, "\n";
    chdir "$path_to_db/$db";
    print "unzipppig ... \n";
    `nice +20 gunzip *.gz`;

    print "making db  ... \n";
    $cmd = "mysqladmin -u root create $db";
    (system $cmd) && warn "error running $cmd\n";

    $cmd = "mysql -u root $db < $db.sql";
    (system $cmd) && warn "error running $cmd\n";

    $cmd = "mysqlimport -u root --fields_escaped_by=\\\\ $db -L *.txt";
    (system $cmd) && die "error running $cmd\n";

    print "\n";	
}
