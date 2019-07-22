#!/usr/bin/perl -w
#

use strict;
use warnings;
my $release_num = 97;

my $path_to_db  = "/storage/databases/ensembl-$release_num/mysql";
(-e $path_to_db) || die "$path_to_db not found\n";
chdir $path_to_db;


my @dbs = split "\n", `ls -d *`;



foreach my $db (@dbs) {

    #next if ( $db =~ /core/);
    # comment this to format compara database (it takes very long time)
    next if ($db =~ /compara/);

    print "************************\n";
    print $db, "\n";
    chdir "$path_to_db/$db";


    print "making db  ... \n";
    my $cmd = "mysqladmin --login-path=ivana  create $db";
    (system $cmd) && warn "error running $cmd\n";

    # accomodate change between mysql versions 5.8 and later in allowed default datetime  values:
    `sed 's/0000-00-00/1000-01-01/g'  $db.sql -i`;

    $cmd = "mysql --login-path=ivana $db < $db.sql";
    (system $cmd) && warn "error running $cmd\n";

    $cmd = "mysqlimport --login-path=ivana --fields_escaped_by=\\\\ $db -L *.txt";
    (system $cmd) && die "error running $cmd\n";

    print "\n";	
}
