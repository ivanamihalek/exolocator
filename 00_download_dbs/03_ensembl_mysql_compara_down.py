#!/usr/bin/env python3
from enum import verify

from config import Config
from download_manager import DownloadManager


class ComparaMysqlDwldManager(DownloadManager):

    def __init__(self, config: Config, verify: bool = True):
        super().__init__(config, verify)
        self.file_type  = "mysql"
        self.local_repo = config.mysql_repo

    def construct_remote_file_level_dir(self, remote_topdir, ensembl_compara_name):
        remote_dir = f"{remote_topdir}/{ensembl_compara_name}"
        return remote_dir

    def construct_local_file_level_dir(self, ensembl_compara_name):
        local_dir = f"{self.local_repo}/{ensembl_compara_name}"
        return local_dir

    def downloadable_files_selection(self, ensembl_compara_name):
        # 'homology_member.txt.gz' has traditionally been the largest, going into 10s of GB
        downloadable = ['homology.txt.gz', f"{ensembl_compara_name}.sql.gz",
                        'ncbi_taxa_name.txt.gz', 'ncbi_taxa_node.txt.gz',
                        'gene_member.txt.gz', 'genome_db.txt.gz',
                        'species_set.txt.gz', 'method_link_species_set.txt.gz',
                        'species_tree_node.txt.gz', 'homology_member.txt.gz']

        return downloadable

    def find_valid_species(self):
        return [f"ensembl_compara_{self.config.release_number}"]


def main():
    #  'homology_member.txt.gz'  can be huge, and thee checksum verification goes on forever
    download_manager = ComparaMysqlDwldManager(Config(), verify=True)
    download_manager.run()


############################################
if __name__ == "__main__":
    main()

