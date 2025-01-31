import functools
import logging
import os
from abc import ABC, abstractmethod
from ftplib import FTP

from config import Config
from download_utils import download_file, download_and_verify, ftp_connect, is_valid_species


class DownloadManager(ABC):

    def __init__(self, config: Config, verify: bool = True):
        self.config = config
        self.ftp_address = "ftp.ensembl.org"
        self.file_type  = None
        self.local_repo = None
        self.ftp: str | FTP = "I am not initialized"
        self.verify_checksum = verify
        # Setup logging
        logging.basicConfig(
            filename='ensembl_download.log',
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s: %(message)s'
        )

    @abstractmethod
    def downloadable_files_selection(self, species):
        pass

    @abstractmethod
    def construct_remote_file_level_dir(self, remote_topdir, species):
        pass

    @abstractmethod
    def construct_local_file_level_dir(self, species):
        pass

    @abstractmethod
    def find_valid_species(self):
        pass

    def download_from_file_level_dir(self, remote_topdir,  species):
        local_dir = self.construct_local_file_level_dir(species)
        os.makedirs(local_dir, exist_ok=True)
        os.chdir(local_dir)

        remote_dir = self.construct_remote_file_level_dir(remote_topdir, species)
        self.ftp.cwd(remote_dir)

        # Download checksums
        checksum_fnm = "CHECKSUMS"
        download_file(self.ftp, checksum_fnm, checksum_fnm)

        # Filter and download files
        downloadable_files = self.downloadable_files_selection(species)

        for item in downloadable_files:
            download_and_verify(self.ftp, item, checksum_fnm, self.verify_checksum)

    ############################
    def run(self):

        if self.file_type is None or self.local_repo is None:
            raise Exception("improperly initialized DownloadManager class")

        # Connect to FTP
        try:
            self.ftp = ftp_connect(self.ftp_address, user='anonymous', passwd='-anonymous@')
        except Exception as e:
            print(e)
            exit()
        # move to our dir
        remote_topdir = f"/pub/release-{self.config.release_number}/{self.file_type}"
        self.ftp.cwd(remote_topdir)

        # Filter valid species using filter() and is_valid_species function
        valid_species = self.find_valid_species()
        for species in valid_species:
            print(f"\nProcessing species: {species}")
            self.download_from_file_level_dir(remote_topdir, species)
            print(f"{species} done.")

        try: self.ftp.quit()
        except: self.ftp.close()
