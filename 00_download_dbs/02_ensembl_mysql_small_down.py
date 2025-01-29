#!/usr/bin/env python3
import functools

from config import Config
from download_manager import DownloadManager
from download_utils import find_unique_file, is_valid_species


class SmallMysqlDwldManager(DownloadManager):

    def __init__(self, config: Config):
        super().__init__(config)
        self.file_type  = "mysql"
        self.local_repo = config.mysql_repo
        self.intermediate_level_dirs = ["mysql"]

    def construct_remote_file_level_dir(self, remote_topdir, species):
        prefix = f"{species}_core_{self.config.release_number}_"
        remote_dir_name = find_unique_file(self.ftp, remote_topdir, prefix)
        remote_dir = f"{remote_topdir}/{remote_dir_name}"
        return remote_dir

    def construct_local_file_level_dir(self, species):
        local_dir = f"{self.local_repo}/{species}_core"
        return local_dir

    def downloadable_files_selection(self, species):
        downloadable = []
        too_big = ['dna.txt.gz', 'repeat_feature.txt.gz', 'CHECKSUMS']
        for file_name in self.ftp.nlst():
            if file_name in too_big: continue
            # CCDS info, contained in dna_align_feature.txt.gz
            # covers confirmed alt splices, but only for human and mouse
            if (file_name in ['dna_align_feature.txt.gz', 'protein_align_feature']
                and species not in ['homo_sapiens', 'mus_musculus']): continue
            downloadable.append(file_name)
        return downloadable

    def find_valid_species(self):
        return list(filter(
            functools.partial(is_valid_species,
                              skip_species=self.config.skip_species,
                              breed_species=self.config.breed_species),
            [dirname.split("_core_")[0] for dirname in self.ftp.nlst() if "_core_" in dirname]
        ))


def main():
    # the intermediate directories here can be "dna", "pep",  or "both"
    download_manager = SmallMysqlDwldManager(Config())
    download_manager.run()


if __name__ == "__main__":
    main()

