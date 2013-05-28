#! /usr/bin/perl -w

@ARGV ||
    die "Usage:  $0  <file name> \n";


@previous_lists = ('actin_binding.txt', 'cell_cycle_checkpoints.txt', 
		   'circadian_clock.txt', 'egfr_signaling.txt', 
		   'enzymes.txt', 'genecards_top500.txt', 
		   'telomere_maintenance.txt', 'wnt_pathway.txt', 
		   'cell_junction.txt', 'transcription.txt', 
		   'translation.txt', 'meiosis.txt');

$filename = $ARGV[0];
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

($blah, $uniprot_id) = ();
%seen = ();

while ( <IF> ) {
    chomp;
    @aux = split "\t";
    for $field(@aux) {
	next if $field !~ /^uniprotkb/;
	($blah, $uniprot_id) = split '\:', $field;
	if ($uniprot_id =~ /\-/) {
	    ($uniprot_id, $blah) = split '\-', $uniprot_id;
	}

	next if $seen{$uniprot_id};

	$seen{$uniprot_id} = 1;
	$ret = `grep $uniprot_id hugo2entrez2uniprot.txt`;

	if (!$ret) {
	    $ensembl_id =  "not_found";
	    next;
	}
	chomp $ret;
	@aux2 = split "\t", $ret;
	
	$ensembl_id = pop @aux2;
	($ensembl_id =~ /^ENSG/) || next;

	$ret = "";
	foreach (@previous_lists) {
	    $ret = `grep $ensembl_id  $_`;
	    last  if $ret;
	}
	next if $ret;

	print "$ensembl_id \n";

	#print "$uniprot_id 	$ensembl_id \n";
    }

}

close IF;
