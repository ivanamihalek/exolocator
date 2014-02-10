#! /usr/bin/perl -w

# I have to be on reindeer and klogged for this to work

use Time::localtime;
use File::stat;

$| = 1; # force stdout flush


$from_dir = "/afs/bii.a-star.edu.sg/dept/biomodel_design/Group/ivana/ExoLocator/results/dumpster/para";
$to_dir   = "/home/ivanam/exolocator/Exolocator/Best_MSA/para";


@species_dirs = split "\n", `ls $from_dir`;

for $spec_dir( @species_dirs ) {

    print  $spec_dir, "\n";

    for $seq_type ( 'pep', 'dna') {

	$target_dir = "$to_dir/$spec_dir/$seq_type";
	(-e $target_dir) || `mkdir -p $target_dir`;

	for $donkey_full_path ( split "\n", `ls $from_dir/$spec_dir/$seq_type/*.afa`) {
=pod
	    @aux = split "/", $donkey_full_path;
	    $afa = pop @aux;
	    $afa || next; # move onif it is an empty string
	    $reindeer_full_path = "$to_dir/$spec_dir/$seq_type/$afa"; 
	    $afa_exists = (  -e  $reindeer_full_path ) ? 1 : 0 ;

	    if ( $afa_exists ) {
		# mtime =   last modify time in seconds since the epoch
		$last_modify_time_on_donkey   = stat($donkey_full_path)->mtime;
		$last_modify_time_on_reindeer = stat($reindeer_full_path)->mtime;

		$but_is_old = ($last_modify_time_on_donkey > $last_modify_time_on_reindeer);
	    }

	    if ( !$afa_exists || $but_is_old) {
		$cmd = "cp $donkey_full_path  $target_dir";
		(system $cmd) && die  "error running $cmd\n";
	    
	    }
=cut
	}
	print "\t $seq_type done\n";
    }
}

