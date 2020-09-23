#! /usr/bin/perl -w

use strict;
use warnings FATAL => 'all';
use Net::FTP;
my $release_num = 101;

#my $local_repository = 
#    "/mnt/ensembl-mirror/release-".$release_num."/mysql";
my $local_repository    = "/storage/databases/ensembl-$release_num/mysql";

(-e  $local_repository) ||
    die "$local_repository not found.\n";

my $ftp_address = "ftp.ensembl.org";

my $ftp = Net::FTP->new( $ftp_address , Debug => 0, Passive=> 1)
    or die "Cannot connect to $ftp_address: $@";

$ftp->login("anonymous",'-anonymous@')
    or die "Cannot login ", $ftp->message;

my $topdir = "/pub/release-$release_num/mysql";

$ftp->cwd($topdir)
    or die "Cannot cwd to $topdir: ", $ftp->message;
$ftp->binary;

my @farm = $ftp->ls;
my $animal;
my ($dir, $local_dir, $foreign_dir,  @contents, $item, $unzipped);

my @skip = ("ancestral_alleles", "caenorhabditis_elegans",
	    "ciona_intestinalis",  "ciona_savignyi", "drosophila_melanogaster",
	    "saccharomyces_cerevisiae");
# ensembl has started collecting breeds for some animals - not of interest here
my @breed = ("mus_musculus_", "hybrid", "oryzias_latipes_", "sus_scrofa_", "capra_hircus_",
			"cyprinus_carpio_", "ovis_aries_");


my @dirs_I_need = ();
my $compara_dir = "ensembl_compara_".$release_num; # take care of it separately

foreach $dir ( @farm ) {
    #$dir =~ /homo_sapiens/ || next;
    ($dir =~ /core/) || next;
    my @aux    = split "_core_", $dir;
    $animal    = $aux[0];

  	next if ( grep {/$animal/} @skip);
	next if ( grep {$animal=~/$_/} @breed);
	# there is no canis_lupus without extension
	next if ($animal=~/canis_lupus_/ &&  $animal ne "canis_lupus_familiaris");
	next if ($animal=~/cricetulus_griseus_/ &&  $animal!~/cricetulus_griseus_crigri/);
	# as of Sept 2020, there are two duck assemblies, both from China:
	# anas_platyrhynchos_platyrhynchos has been updated more recently
	# and has more gene transcripts then anas_platyrhynchos
	next if ($animal=~/anas_platyrhynchos/ &&  $animal!~/anas_platyrhynchos_platyrhynchos/);

    push @dirs_I_need, $dir;
}


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
	
    @contents =  $ftp->ls;

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
            print LOG   "getting $item  failed ", $ftp->message, "\n";
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

#
#################################################################################
# now take care of compara db
# caveat here: homology_member.txt.gz is huge - almo 70 GB as og ensembl 101
# 15 hours to download (using wget), another 1 hr to gunzip to 173 GB
# can I reconstruct the contents on my own, without downloading?
# the time to load this shit into mySQL is also becoming ridiculous  - 4hr just for homology table
$dir       = $compara_dir;
$local_dir = "$local_repository/$dir" ;
(-e $local_dir )|| `mkdir -p $local_dir`;

$foreign_dir = "$topdir/$dir";

$ftp->cwd($foreign_dir)
    or die "Cannot cwd to $foreign_dir: ", $ftp->message;
@contents =  $ftp->ls;

my $compara_sql = 'ensembl_compara_'.$release_num.'.sql.gz';
 
foreach $item ('homology.txt.gz', 'homology_member.txt.gz', 
	       'gene_member.txt.gz', 'genome_db.txt.gz',
            'species_set.txt.gz', 'method_link_species_set.txt.gz',
            'species_tree_node.txt.gz', $compara_sql) {

    (grep {/$item/} @contents) || die "$item not found in $foreign_dir\n";
    
    print "\t$item\n";
    $unzipped = $item;
    $unzipped =~ s/\.gz$//;
    if ( -e "$local_dir/$unzipped" ) {
    	print "\t\t $unzipped found in $local_dir\n";
    	next;
    }

    if ( ! $ftp->get($item) ) {
    	print LOG   "getting $item  failed ", $ftp->message, "\n";
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
# Paired-end tags (PET) (sometimes "Paired-End diTags", or simply "ditags") are the short sequences
#at the 5’ and 3’ ends of a DNA fragment which are unique enough that they (theoretically)
#exist together only once in a genome, therefore making the sequence of the DNA in between them available upon search
#  - we are not using that feature (should or could we?)
