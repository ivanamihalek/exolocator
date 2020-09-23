#! /usr/bin/perl -w

# TODO: for the next update - download species names
# and decide which of them should be downloaded prior to downloading
# - check if they correspond to different taxonomy_ids, and which build is nerwe
# apparently ensembl started shoveling
# uncurated number of builds so we have
# astyanax_mexicanus_pachon is an older build of astyanax_mexicanus
# and anas_platyrhynchos and anas_platyrhynchos_platyrhynchos
# which are not subspecies but different builds
# I am not sure what use if any I have from haveing both male and female H.glaber
# male annotation is newer

# notes:
# *pep.all.fa vs *pep.abinitio.fa
# not clear if there is any overlap. Genes are labeled with GENSCAN rather than ENS.
# THese files, however are 10 to 20 times smaller than dna files, so it is not clear
# if there is mauch gain in not downloading them/

# Dna seqeunces: go for rm toplevel versions:
# rm stands for repat masked, and toplevel stands for "all sequence regions flagged as toplevel in an Ensembl
# schema. This includes chromsomes, regions not assembled into chromosomes and
# N padded haplotype/patch regions."

# TODO: check for the files that exist locally in the older release directory
# (or should I just delete everything except human once I am done putting together the ?)

use strict;
use warnings FATAL => 'all';
use Net::FTP;
my $release_num = 101;
my $local_repository = "/storage/databases/ensembl-$release_num/fasta";
$| = 1; # flush stdout

-e $local_repository ||
    die "local repository:\n$local_repository\nnot found\n";

my $ftp_address = "ftp.ensembl.org";

my $ftp = Net::FTP->new( $ftp_address , Debug => 0, Passive=> 1)
    or die "Cannot connect to $ftp_address: $@";

$ftp->login("anonymous",'-anonymous@')
    or die "Cannot login ", $ftp->message;

my $topdir = "/pub/release-$release_num/fasta";
$ftp->cwd($topdir)
    or die "Cannot cwd to $topdir: ", $ftp->message;
$ftp->binary;

my @farm = $ftp->ls;
my $animal;

my @skip = ("ancestral_alleles", "caenorhabditis_elegans",
	    "ciona_intestinalis",  "ciona_savignyi", "drosophila_melanogaster",
	    "saccharomyces_cerevisiae");
# ensembl has started collecting breeds for some animals - not of interest here
my @breed = ("mus_musculus_", "hybrid", "oryzias_latipes_", "sus_scrofa_", "capra_hircus_",
			"cyprinus_carpio_", "ovis_aries_", "astyanax_mexicanus_");
# astyanax_mexicanus_pachon is an older build of astyanax_mexicanus

my ($dir, $local_dir, $foreign_dir,  @contents, $item, $unzipped);

open (LOG, ">ensembl_download.log") || die "error opening log: $!\n";

my $ct = 0;

#@farm = ('homo_sapiens');
foreach $animal ( @farm ) {

	#print("$animal\n");
	next if ( grep {/$animal/} @skip);
	next if ( grep {$animal=~/$_/} @breed);
	# there is no canis_lupus without extension
	next if ($animal=~/canis_lupus_/ &&  $animal ne "canis_lupus_familiaris");
	next if ($animal=~/cricetulus_griseus_/ &&  $animal!~/cricetulus_griseus_crigri/);
	# as of Sept 2020, there are two duck assemblies, both from China:
	# anas_platyrhynchos_platyrhynchos has been updated more recently
	# and has more gene transcripts then anas_platyrhynchos
	next if ($animal=~/anas_platyrhynchos/ &&  $animal!~/anas_platyrhynchos_platyrhynchos/);

    $ct += 1;
    print $ct, "  ", $animal, "\n";


    foreach $dir  ( "pep",  "dna" ){

		$local_dir = "$local_repository/$animal/$dir" ;
		( -e $local_dir )  || `mkdir -p $local_dir`;

		$foreign_dir = "$topdir/$animal/$dir";

		$ftp->cwd($foreign_dir)
		    or die "Cannot cwd to $foreign_dir: ", $ftp->message;

		# download checksums
		my $checksums= "CHECKSUMS";
		$ftp->get($checksums) || die "getting $checksums  failed:", $ftp->message;

		@contents =  $ftp->ls;
		foreach $item (@contents) {

		    next if ($item!~/\.dna_rm\.toplevel.*gz$/ && $item!~/.*pep.*gz$/ );

		    print LOG "\t$item\n";
			print "\t$item\n";

		    $unzipped = $item;
		    $unzipped =~ s/\.gz$//;
		    if ( -e "$local_dir/$unzipped" ) {
				print  LOG  "\t\t $unzipped found in $local_dir\n";
				next;
		    }

		    $ftp->get($item) || die "getting $item  failed ", $ftp->message;
			# do the checsum
			my $sum_reported = `grep $item $checksums | awk '{print \$1, \$2}'`; chomp $sum_reported;
			print "\t\t$item downloaded, running checksum\n";
			my $sum = `sum $item | awk '{print \$1, \$2}'`; chomp $sum;
			($sum_reported eq $sum) || die " checksum mismatch\n";

		    `mv  $item  $local_dir`;
	    
		    print  LOG "\t\t $item moved to $local_dir\n";

		    if (system ( "gunzip $local_dir/$item" )) {
				print LOG "error uncompressing $local_dir/$item.\n";
				die  "\t\terror uncompressing $local_dir/$item.\n";
		    } else {
				print "\t\t $item unzipped \n";
		    }
		}

		# remove the old checksums file
		`rm -f $checksums`;
    }
}
   
close LOG;
