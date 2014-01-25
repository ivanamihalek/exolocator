#! /usr/bin/perl -w

use Time::localtime;
use File::stat;

$from_dir = "/afs/bii.a-star.edu.sg/dept/biomodel_design/Group/ivana/ExoLocator/results/dumpster";
$to_dir   = "/home/ivanam/exolocator/Exolocator/Best_MSA";

$filename = "";
@ARGV && ($filename=$ARGV[0]);

if (!$filename or $filename eq 'none') {    
    $home = `pwd`; chomp $home;
    chdir "$from_dir/pep/";
    @ens_ids = split "\n", `ls *afa | sed  's/\.afa//g'`;
    chdir $home; 
} else {
    $filename = $ARGV[0];
    @ens_ids  = split "\n", `awk '{print \$1}' $filename`;
}


foreach  $ensembl_id (@ens_ids) {

    print  $ensembl_id, "\n";

    $donkey_full_path   = "$from_dir/pep/$ensembl_id.afa";
    $reindeer_full_path = "$to_dir/pep/$ensembl_id.afa";

    ( -z $donkey_full_path ) && next;

    $afa_exists    = (  -e  $reindeer_full_path ) ? 1 : 0 ;
    $is_up_to_date = 0;
    if ( $afa_exists ) {
	# mtime =   last modify time in seconds since the epoch
	$last_modify_time_on_donkey   = stat($donkey_full_path)->mtime;
	$last_modify_time_on_reindeer = stat($reindeer_full_path)->mtime;

	$is_up_to_date = ($last_modify_time_on_donkey  <=  $last_modify_time_on_reindeer);
    }

    $afa_exists &&  $is_up_to_date  && next;


    $cmd = "cp $from_dir/pep/$ensembl_id.afa $to_dir/pep/$ensembl_id.afa";
    print $cmd, "\n";
    (system $cmd) && print "error ...\n";
	
    $cmd = "cp $from_dir/dna/$ensembl_id.afa $to_dir/dna/$ensembl_id.afa";
    print $cmd, "\n";
    (system $cmd) && print "error ...\n";
	
    $cmd = "cp $from_dir/notes/$ensembl_id.txt $to_dir/notes/$ensembl_id.txt";
    print $cmd, "\n";
    (system $cmd) && print "error ...\n";

    print `grep SW  $to_dir/notes/$ensembl_id.txt | wc -l`;
    print "\n";
	
}

