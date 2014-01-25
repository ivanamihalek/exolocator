#! /usr/bin/perl -w

@ARGV ||
    die "Usage:  $0  <file name> \n";

$from_dir = "/afs/bii.a-star.edu.sg/dept/biomodel_design/Group/ivana/ExoLocator/results/dumpster";
$to_dir   = "/home/ivanam/exolocator/Exolocator/Best_MSA";


if (!$filename or $filename eq 'none') {
    $home = `pwd`; chomp $home;
    chdir "$from_dir/pep/";
    @ens_ids = split "\n", `ls *afa | sed  's/\.afa//g'`;
    chdir $pwd; 
} else {
    $filename = $ARGV[0];
    @ens_ids = split "\n", `awk '{print \$1}' $filename`;
}

print "@ens_ids[0..9]\n";
exit(1);

foreach  $ensembl_id ( @ens_ids) {
    print  $ensembl_id, "\n";

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

close IF;
