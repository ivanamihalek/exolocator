#!/usr/bin/env python3

from config import Config
from download_manager import DownloadManager


class FastaDwldManager(DownloadManager):
    
    def __init__(self, config: Config):
        super().__init__(config)
        self.file_type  = "fasta"
        self.local_repo = config.fasta_repo
        # the intermediate directories here can be "dna", "pep",  or "both"
        self.intermediate_level_dirs = ["pep"]

    def downloadable_files_selection(self, species):
        # species is not really needed here, we needed it when implemetnign this function for mysql files
        return list(filter(
            lambda x: x.endswith(('.dna.toplevel.gz', '.pep.all.fa.gz')), self.ftp.nlst()
        ))


def main():

    download_manager = FastaDwldManager(Config())
    download_manager.run()


if __name__ == "__main__":
    main()
