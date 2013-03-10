#!/usr/bin/perl -w
# one extra database that is needed,
# never made it a part of this script:
# ensembl_compara_xx, the tables needed are homology, homology_member, member and genome_db
# the three largest chunks in *core*: dna, dna_align_feature, and protein_alig_feature
# (actually, from dna_align_feature can be found the info about CCDS annotation)
# I nevere really found use for - perhaps can be skipped during downloading, to make things faster
# the sql script needs to be hacked (or does it - just leave the tables empty)

$path_to_db = "/mnt/ensembl/release-69/mysql";

$credentials = " -h jupiter.private.bii -P 3307 -u root -psqljupitersql  ";

(-e $path_to_db) || die "$path_to_db not found\n";

chdir $path_to_db;


@dbs = split "\n", `ls -d *`;



foreach $db (@dbs) {

    next if ( $db =~ /core/); # <<<<< negate this to format compara database


    print "************************\n";
    print $db, "\n";
    chdir "$path_to_db/$db";
    print "unzipppig ... \n";
    `nice +20 gunzip *.gz`;

    print "making db  ... \n";
    $cmd = "mysqladmin $credentials  create $db";
    (system $cmd) && warn "error running $cmd\n";

    $cmd = "mysql $credentials $db < $db.sql";
    (system $cmd) && warn "error running $cmd\n";

    $cmd = "mysqlimport $credentials --fields_escaped_by=\\\\ $db -L *.txt";
    (system $cmd) && die "error running $cmd\n";

    print "\n";	
}