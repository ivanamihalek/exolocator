#! /usr/bin/perl -w

@ARGV ||
    die "Usage:  $0  <file name> \n";

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

	if ($ret) {
	    chomp $ret;
	    @aux2 = split "\t", $ret;
	
	    $ensembl_id = pop @aux2;
	    ($ensembl_id =~ /^ENSG/) && print "$ensembl_id \n";

	} else {

	    $ensembl_id =  "not_found";
	}
	#print "$uniprot_id 	$ensembl_id \n";
    }

}

close IF;
