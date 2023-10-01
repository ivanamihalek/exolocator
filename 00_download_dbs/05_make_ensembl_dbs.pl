#!/usr/bin/perl -w
#  useful here:
#  https://www.digitalocean.com/community/tutorials/how-to-move-a-mysql-data-directory-to-a-new-location-on-ubuntu-18-04

use strict;
use warnings;
my $release_num = 110;

my $path_to_db  = "/storage/databases/ensembl-$release_num/mysql";
(-e $path_to_db) || die "$path_to_db not found\n";
chdir $path_to_db;

my $dry = 0;

my @dbs = split "\n", `ls -d *`;

for (my $idx = 0; $idx < scalar(@dbs); $idx++) {
    my $db = $dbs[$idx];

    #    next if ( $db =~ /^[abc]/);
    #    next if ( $db eq "danio_rerio_core_101_11");
    # comment this to format compara database (it takes very long time)
    # next if ($db !~ /compara/);

    print "************************\n";
    print $db, "\n";
    chdir "$path_to_db/$db";


    # setting longpath: mysql_config_ed
    # https://dev.mysql.com/doc/refman/5.6/en/mysql-config-editor.html
    # in mysql shell: SET GLOBAL show_compatibility_56 = ON;
    # see: http://mysql-tools.com/home/1-mysql/125-showcompatibility56-mysql-57-connection-problem.html
    print "making db  ... \n";
    my $cmd = "mysqladmin --login-path=tcga  create $db";
    if ($dry) {
        print($cmd, "\n");
    } else {
        (system $cmd) && die "error running $cmd\n";
    }

    ########################################
    (-e "$db.sql.gz" ) || die "$db.sql.gz not found\n";
    `pigz -d -k $db.sql.gz`;  # k keeps the original file
    # accommodate change between mysql versions 5.8 and later in allowed default datetime  values:
    `sed 's/0000-00-00/1000-01-01/g'  $db.sql -i`;

    # as of mysql v 8,  NFORMATION_SCHEMA.SESSION_VARIABLES became performance_schema.session_variables
    `sed 's/INFORMATION_SCHEMA.SESSION_VARIABLES/performance_schema.session_variables/g' -i  $db.sql `;

    # TODO load without gunzipping
    # zcat /path/to/file.sql.gz | mysql -u 'root' -p your_database
    $cmd = "mysql --login-path=tcga $db < $db.sql";
    if ($dry) {
        print($cmd, "\n");
    } else {
        (system $cmd) && warn "error running $cmd\n";
    }
    `rm $db.sql`;


    my @txt_files = ();
    for my $txtf ( glob("*.txt.gz")) {
        `pigz -d -k $txtf`;
        $txtf =~ s/\.gz$//;
        push @txt_files, $txtf
    }

    $cmd = "mysqlimport --login-path=tcga --fields_escaped_by=\\\\ $db -L *.txt";
    if ($dry) {
        print($cmd, "\n");
    } else {
        (system $cmd) && die "error running $cmd\n";
    }

     for my $txtf (@txt_files) {
        `rm $txtf`;
     }

    print "\n";
    exit;
}
