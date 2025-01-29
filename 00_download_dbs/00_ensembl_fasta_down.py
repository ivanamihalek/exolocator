#!/usr/bin/env python3
import functools
import os

from config import Config
from download_manager import DownloadManager
from download_utils import is_valid_species


class FastaDwldManager(DownloadManager):
    
    def __init__(self, config: Config, fasta_file_type):
        super().__init__(config)
        self.file_type  = "fasta"
        self.local_repo = config.fasta_repo
        # the fast_file_type here can be "dna", "pep",  or "both"
        self.fasta_file_type = fasta_file_type

    def construct_remote_file_level_dir(self, remote_topdir, species):
        remote_dir = f"{remote_topdir}/{species}/{self.fasta_file_type}"
        return remote_dir

    def construct_local_file_level_dir(self, species):
        local_dir = os.path.join(self.local_repo, species, self.fasta_file_type)
        return local_dir

    def downloadable_files_selection(self, species):
        # species is not really needed here, we needed it when implemetnign this function for mysql files
        return list(filter(
            lambda x: x.endswith(('.dna.toplevel.gz', '.pep.all.fa.gz')), self.ftp.nlst()
        ))

    def find_valid_species(self):
        return list(filter(
            functools.partial(is_valid_species,
                              skip_species=self.config.skip_species,
                              breed_species=self.config.breed_species),
            self.ftp.nlst()
        ))


def main():
    # to do the actual exon patching, we'll need both "dna: and "pep" fasta files
    # to construct alignments of canonical peptides, pep should be enough
    download_manager = FastaDwldManager(Config(), "pep")
    download_manager.run()


if __name__ == "__main__":
    main()
