#!/usr/bin/env python3

from config import Config
from download_manager import DownloadManager


class SmallMysqlDwldManager(DownloadManager):

    def __init__(self, config: Config):
        super().__init__(config)
        self.file_type  = "mysql"
        self.local_repo = config.mysql_repo
        self.intermediate_level_dirs = ["mysql"]

    def downloadable_files_selection(self, species):
        downloadable = []
        too_big = ['dna.txt.gz', 'repeat_feature.txt.gz']
        for file_name in self.ftp.nlst():
            if file_name in too_big: continue
            # CCDS info, contained in dna_align_feature.txt.gz
            # covers confirmed alt splices, but only for human and mouse
            if (file_name in ['dna_align_feature.txt.gz', 'protein_align_feature']
                and species not in ['homo_sapiens', 'mus_musculus']): continue
            downloadable.append(file_name)
        return downloadable


def main():
    # the intermediate directories here can be "dna", "pep",  or "both"
    download_manager = SmallMysqlDwldManager(Config())
    download_manager.run()


if __name__ == "__main__":
    main()

