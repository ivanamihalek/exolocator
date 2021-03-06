#! /usr/bin/perl -w


@ARGV==3 || die "Usage: $0 <set name> <number of threads> <method>.\n";

#($set, $no_thr, $method) = @ARGV;



foreach (  "18_find_missing_exons.py",    "20_novel_exon_cleanup.py",
	   "17_ortho_exon_map_to_msa.py", "07_gene2exon_store.py",
	   "21_make_novel_exon_maps.py",  "17_ortho_exon_map_to_msa.py",
	   "30_db_migration/06_reconstruct_ortho_alnmts.py") {

    -e $_ || die "$_ not found.\n";
}
($set, $no_threads, $method) =  @ARGV;

#18_find_missing_exons
$cmd = "18_find_missing_exons.py  $set  $no_threads  $method ";
(system $cmd) && die "error running $cmd";


#20_novel_exon_cleanup
$cmd = "20_novel_exon_cleanup.py   $set  $no_threads ";
(system $cmd); # && die "error running $cmd";

#07_gene2exon_store 
$cmd = "07_gene2exon_store.py   $set  $no_threads ";
(system $cmd) && die "error running $cmd";


#21_make_novel_exon_maps
$cmd = "21_make_novel_exon_maps.py   $set  $no_threads ";
(system $cmd) && die "error running $cmd";

#17_ortho_exon_map_to_msa
$cmd = "17_ortho_exon_map_to_msa.py   $set  $no_threads ";
(system $cmd) && die "error running $cmd";


#30/06     alignement reconstruction
$cmd = "30_db_migration/06_reconstruct_ortho_alnmts.py   $set  $no_threads ";
(system $cmd) && die "error running $cmd";
