#! /usr/bin/perl -w

# notes:
# *pep.all.fa vs *pep.abinitio.fa
# not clear if there is any overlap. Genes are labeled with GENSCAN rather than ENS.
# THese files, however are 10 to 20 times smaller than dna files, so it is not clear
# if there is mauch gain in not downloading them/

# Dna seqeunces: go for rm toplevel versions:
# rm stands for repat masked, and toplevel stands for "all sequence regions flagged as toplevel in an Ensembl
# schema. This includes chromsomes, regions not assembled into chromosomes and
# N padded haplotype/patch regions."



use strict;
use warnings FATAL => 'all';
use Net::FTP;
my $release_num = 97;
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

my ($dir, $local_dir, $foreign_dir,  @contents, $item, $unzipped);

open (LOG, ">ensembl_download.log") || die "error opening log: $!\n";

my $ct = 0;

#@farm = ('homo_sapiens');
foreach $animal ( @farm ) {

	#print("$animal\n");
    next if ( grep {/$animal/} @skip);
  	# there is a bunch of mouse genomes that I do not want to deal with now
	next if ($animal=~/mus_musculus_/ || $animal=~/hybrid/ || $animal=~/oryzias_latipes_/ || $animal=~/sus_scrofa_/);

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
