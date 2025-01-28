import functools
import logging
import os
from abc import ABC, abstractmethod
from ftplib import FTP

from config import Config
from download_utils import download_file, download_and_verify, ftp_connect, is_valid_species


class DownloadManager(ABC):

    def __init__(self, config: Config):
        self.config = config
        self.ftp_address = "ftp.ensembl.org"
        self.file_type  = None
        self.local_repo = None
        self.intermediate_level_dirs = []
        self.ftp: str | FTP = "I am not intialized"
        # Setup logging
        logging.basicConfig(
            filename='ensembl_download.log',
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s: %(message)s'
        )

    @abstractmethod
    def downloadable_files_selection(self, species):
        pass

    def download_from_file_level_dir(self, species, foreign_topdir, file_level_dir):
        local_dir = os.path.join(self.local_repo, species, file_level_dir)
        os.makedirs(local_dir, exist_ok=True)
        os.chdir(local_dir)

        foreign_dir = f"{foreign_topdir}/{species}/{file_level_dir}"
        self.ftp.cwd(foreign_dir)

        # Download checksums
        checksum_fnm = "CHECKSUMS"
        download_file(self.ftp, checksum_fnm, checksum_fnm)

        # Filter and download files
        downloadable_files = self.downloadable_files_selection(species)

        for item in downloadable_files:
            download_and_verify(self.ftp, item, checksum_fnm)

    def species_download(self,  animal, foreign_topdir):
        for file_level_dir in self.intermediate_level_dirs:
            self.download_from_file_level_dir(animal, foreign_topdir, file_level_dir)

    ############################
    def run(self):

        if self.file_type is None or self.local_repo is None or  not self.intermediate_level_dirs:
            raise Exception("improperly initialized DownloadManager class")

        # Connect to FTP
        try:
            self.ftp = ftp_connect(self.ftp_address, user='anonymous', passwd='-anonymous@')
        except Exception as e:
            print(e)
            exit()
        # move to our dir
        foreign_topdir = f"/pub/release-{self.config.release_number}/{self.file_type}"
        self.ftp.cwd(foreign_topdir)

        # Filter valid species using filter() and is_valid_species function
        valid_species = list(filter(
            functools.partial(is_valid_species,
                              skip_species=self.config.skip_species,
                              breed_species=self.config.breed_species),
            self.ftp.nlst()
        ))

        for species in valid_species:
            print(f"\nProcessing species: {species}")
            self.species_download(species, foreign_topdir)
            print(f"{species} done.")

        try: self.ftp.quit()
        except: self.ftp.close()
