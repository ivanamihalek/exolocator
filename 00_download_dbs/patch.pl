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


my ($dir, $local_dir, $foreign_dir,  @contents, $item, $unzipped);


open (LOG, ">ensembl_mysql_download.log") || die "error opening log: $!\n";

#################################################################################
my $ct = 0;

foreach $dir ( 'homo_sapiens_core_101_38' , 'mus_musculus_core_101_38')  {

    $ct += 1;
    print "$ct   $dir \n";

    $local_dir = "$local_repository/$dir" ;
    (-e $local_dir )|| `mkdir -p $local_dir`;

    $foreign_dir = "$topdir/$dir";

    $ftp->cwd($foreign_dir)
	or die "Cannot cwd to $foreign_dir: ", $ftp->message;
	
    @contents =  ('protein_align_feature.txt.gz');

    foreach $item (@contents) {

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

close LOG;



###################################
# Paired-end tags (PET) (sometimes "Paired-End diTags", or simply "ditags") are the short sequences
#at the 5’ and 3’ ends of a DNA fragment which are unique enough that they (theoretically)
#exist together only once in a genome, therefore making the sequence of the DNA in between them available upon search
#  - we are not using that feature (should or could we?)
