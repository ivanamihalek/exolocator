#!/usr/bin/env python3
import functools

from config import Config
from download_manager import DownloadManager
from download_utils import find_unique_file, is_valid_species


class ComparaMysqlDwldManager(DownloadManager):

    def __init__(self, config: Config):
        super().__init__(config)
        self.file_type  = "mysql"
        self.local_repo = config.mysql_repo
        self.intermediate_level_dirs = ["mysql"]

    def construct_remote_file_level_dir(self, remote_topdir, ensembl_compara_name):
        remote_dir = f"{remote_topdir}/{ensembl_compara_name}"
        return remote_dir

    def construct_local_file_level_dir(self, ensembl_compara_name):
        local_dir = f"{self.local_repo}/{ensembl_compara_name}"
        return local_dir

    def downloadable_files_selection(self, ensembl_compara_name):
        downloadable = ['homology.txt.gz', 'homology_member.txt.gz',
                        'gene_member.txt.gz', 'genome_db.txt.gz',
                        'species_set.txt.gz', 'method_link_species_set.txt.gz',
                        'species_tree_node.txt.gz', f"{ensembl_compara_name}.sql.gz"]

        return downloadable

    def find_valid_species(self):
        return [f"ensembl_compara_{self.config.release_number}"]


def main():
    download_manager = ComparaMysqlDwldManager(Config())
    download_manager.run()


if __name__ == "__main__":
    main()

