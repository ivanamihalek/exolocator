#! /usr/bin/perl -w

use strict;
use Net::FTP;
my $release_num = 86;

#my $local_repository = 
#    "/mnt/ensembl-mirror/release-".$release_num."/mysql";
my $local_repository    = "/databases/ensembl-$release_num/mysql";

(-e  $local_repository) ||
    die "$local_repository not found.\n";

my $ftp_address = "ftp.ensembl.org";

my $ftp = Net::FTP->new( $ftp_address , Debug => 0, Passive=> 1)
    or die "Cannot connect to $ftp_address: $@";

$ftp->login("anonymous",'-anonymous@')
    or die "Cannot login ", $ftp->message;

my $topdir = "/pub/release-".$release_num."/mysql";

$ftp->cwd($topdir)
    or die "Cannot cwd to $topdir: ", $ftp->message;
$ftp->binary;

my @farm = $ftp->ls;
my $animal;
my ($dir, $local_dir, $foreign_dir,  @contents, $item, $unzipped);


my @skip = ("ancestral_alleles", "caenorhabditis_elegans", "ciona_intestinalis",  
	    "ciona_savignyi", "drosophila_melanogaster", "saccharomyces_cerevisiae");

my @dirs_I_need = ();
my $compara_dir = "ensembl_compara_".$release_num; # take care of it separately

foreach $dir ( @farm ) {
    #$dir =~ /homo_sapiens/ || next;
    ($dir =~ /core/) || next;
    ($dir =~ /expression/) && next;
    my @aux    = split "_", $dir;
    $animal    = join  "_", @aux[0..1];
    (grep {/$animal/} @skip)  && next;
    next if ($animal=~/mus_musculus_/);
    push @dirs_I_need, $dir;
}

print join (" ", @dirs_I_need), "\n";


open (LOG, ">ensembl_mysql_download.log") || die "error opening log: $!\n";

#################################################################################
my $ct = 0;

foreach $dir ( @dirs_I_need) {

    $ct += 1;
    print "$ct   $dir \n";

    $local_dir = "$local_repository/$dir" ;
    (-e $local_dir )|| `mkdir -p $local_dir`;

    $foreign_dir = "$topdir/$dir";

    $ftp->cwd($foreign_dir)
	or die "Cannot cwd to $foreign_dir: ", $ftp->message;
	
    my @contents =  $ftp->ls;
    my $item;

    foreach $item (@contents) {

	next if ($item !~ /\.gz$/);
	next if ($item eq 'assembly.txt.gz');
	next if ($item eq 'dna.txt.gz');
	next if ($item eq 'protein_align_feature.txt.gz');
	next if ($item eq 'repeat_feature.txt.gz');
	next if ($item eq 'ditag_feature.txt.gz'); # see below for ditag definition
	next if ($item eq 'ditag.txt.gz');
	# CCDS info, contained in dna_align_feature.txt.gz  covers confirmed alt splices,
	# but only for human and mouse
	if ($item eq 'dna_align_feature.txt.gz') {
	    next if ($dir !~ 'homo_sapiens' && $dir !~ 'mus_musculus');
	}
	print "\t$item\n";

	$unzipped = $item;
	$unzipped =~ s/\.gz$//;
	if ( -e "$local_dir/$unzipped" ) {
	    print "\t\t $unzipped found in $local_dir\n";
	    next;
	}

	if ( ! $ftp->get($item) ) {
	    print LOG   "getting $item  failed ", $ftp->message;
	    next;
	}

	`mv  $item  $local_dir`;
	    
	print "\t\t $item moved to $local_dir\n";

	if (system ( "gunzip $local_dir/$item" )) {
	    print LOG "error uncompressing $local_dir/$item.\n";
	    print     "\t\terror uncompressing $local_dir/$item.\n";
	} else {
	    print "\t\t $item unzipped \n";
	}
    }
}


=pod
#################################################################################
# now take care of compara db
$dir       = $compara_dir;
$local_dir = "$local_repository/$dir" ;
(-e $local_dir )|| `mkdir -p $local_dir`;

$foreign_dir = "$topdir/$dir";

$ftp->cwd($foreign_dir)
    or die "Cannot cwd to $foreign_dir: ", $ftp->message;
@contents =  $ftp->ls;

my $compara_sql = 'ensembl_compara_'.$release_num.'.sql.gz';
 
foreach $item ('homology.txt.gz', 'homology_member.txt.gz', 
	       'gene_member.txt.gz', 'genome_db.txt.gz', $compara_sql) {

    (grep {/$item/} @contents) || die "$item not found in $foreign_dir\n";
    
    print "\t$item\n";
    $unzipped = $item;
    $unzipped =~ s/\.gz$//;
    if ( -e "$local_dir/$unzipped" ) {
	print "\t\t $unzipped found in $local_dir\n";
	next;
    }

    if ( ! $ftp->get($item) ) {
	print LOG   "getting $item  failed ", $ftp->message;
	next;
    }

    `mv  $item  $local_dir`;
	    
    print "\t\t $item moved to $local_dir\n";

    if (system ( "gunzip $local_dir/$item" )) {
	print LOG "error uncompressing $local_dir/$item.\n";
	print     "\t\terror uncompressing $local_dir/$item.\n";
    } else {
	print "\t\t $item unzipped \n";
    }


}

   
close LOG;



###################################
# Paired-end tags (PET) (sometimes "Paired-End diTags", or simply "ditags") are the short sequences at the 5’ and 3’ ends of a DNA fragment which are unique enough that they (theoretically) exist together only once in a genome, therefore making the sequence of the DNA in between them available upon search
#  - we are not using that feature (should or could we?)
