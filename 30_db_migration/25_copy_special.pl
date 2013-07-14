#! /usr/bin/perl -w

@ARGV ||
    die "Usage:  $0  <file name> \n";

$from_dir = "/afs/bii.a-star.edu.sg/dept/biomodel_design/Group/ivana/ExoLocator/results/dumpster";
$to_dir   = "/home/ivanam/exolocator/Exolocator/Best_MSA";


$filename = $ARGV[0];
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";


while ( <IF> ) {
    chomp;
    @aux = split;
    $ensembl_id = $aux[0];
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
