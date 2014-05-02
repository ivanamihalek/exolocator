#! /usr/bin/perl -w


@ARGV==3 || die "Usage: $0 <set name> <number of threads> <method>.\n";

#($set, $no_thr, $method) = @ARGV;

@scripts = ("09_make_ortho_maps.py",
	    "10_ortho_exon_map_to_msa.py",
	    "11_find_missing_exons.py",
	    "12_novel_exon_cleanup.py",
	    "13_gene2exon_store.py",
	    "14_make_novel_exon_maps.py",
	    "15_ortho_exon_map_to_msa.py",
	    "30_db_migration/06_reconstruct_ortho_alnmts.py");

foreach (@scripts) {

    -e $_ || die "$_ not found.\n";
}
($set, $no_threads, $method) =  @ARGV;

foreach (@scripts) {
    $cmd = "$_  $set  $no_threads  $method";
    (system $cmd) && die "error running $cmd\n";
}
