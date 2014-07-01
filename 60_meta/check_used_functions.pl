#! /usr/bin/perl -w
# change only the scripts that are going to be used immediately


@ARGV ||
    die "Usage:  $0  <file name> \n";

$filename = $ARGV[0];
open (IF, "<$filename" ) 
    || die "Cno $filename: $!.\n";

while ( <IF> ) {
    if (/^def/ ) {
	chomp;
	$line    =  $_;
	$line    =~ s/[\(\)]/ /g;
	@aux     =  split ' ', $line;
	$fn_name =  $aux[1];
	next if $fn_name eq 'main';

	# is this function referenced in this file?
	print "\n===================\n";
	print $fn_name, "\n";
	@retlines = split "\n", `grep -n $fn_name $filename`;
	for $retline (@retlines) {
	    next if $retline =~ /def /;
	    print "** $retline\n";
	}

	# do we have another (or perhaps the same?)
	# function of the same name defined some place else?
	
    }

}

close IF;
