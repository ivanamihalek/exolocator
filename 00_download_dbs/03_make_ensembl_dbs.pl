#!/usr/bin/perl -w
# one extra database that is needed,
# never made it a part of this script:
# ensembl_compara_xx, the tables needed are homology, homology_member, member and genome_db
# the three largest chunks in *core*: dna, dna_align_feature, and protein_alig_feature
# (actually, from dna_align_feature can be found the info about CCDS annotation)
# I nevere really found use for - perhaps can be skipped during downloading, to make things faster
# the sql script needs to be hacked (or does it - just leave the tables empty)


$release_no = 86;

$path_to_db  = "/databases/ensembl-$release_no/mysql";
$credentials = "";

(-e $path_to_db) || die "$path_to_db not found\n";

chdir $path_to_db;


@dbs = split "\n", `ls -d *`;



foreach $db (@dbs) {

    #next if ( $db =~ /core/); # <<<<< uncomment this to format compara database (it takes very long time)
    next if ($db =~ /compara/);

    print "************************\n";
    print $db, "\n";
    chdir "$path_to_db/$db";


    print "making db  ... \n";
    $cmd = "mysqladmin --login-path=client  create $db";
    (system $cmd) && warn "error running $cmd\n";

    # accomodate change between mysql versions 5.8 and later in allowed default datetime  values:
    `sed 's/0000-00-00/1000-01-01/g'  $db.sql -i`;

    $cmd = "mysql --login-path=client $db < $db.sql";
    (system $cmd) && warn "error running $cmd\n";

    $cmd = "mysqlimport --login-path=client --fields_escaped_by=\\\\ $db -L *.txt";
    (system $cmd) && die "error running $cmd\n";

    print "\n";	
}
