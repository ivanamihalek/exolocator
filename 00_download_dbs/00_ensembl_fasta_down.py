#!/usr/bin/env python3

import os
import ftplib
import logging
import functools
import subprocess

from config import Config

def ftp_connect(host, user, passwd):
    ftp = ftplib.FTP(host)
    try:
        ftp.login(user, passwd)
        # FTP operations
    except ftplib.all_errors as ftp_error:
        logging.error(f"FTP Operation Error: {ftp_error}")
        exit(1)
    print(f"logged into {host} as {user}.")
    return ftp


def calculate_checksum(filepath):
    """
    Calculate checksum using Linux shell 'sum' command

    Args:
        filepath (str): Path to the file

    Returns:
        str: Checksum and block count from 'sum' command
    """
    try:
        result = subprocess.run(
            ['sum', filepath],
            capture_output=True,
            text=True,
            check=True
        )
        return result.stdout.split()[0]
    except subprocess.CalledProcessError as e:
        raise ValueError(f"Checksum calculation failed: {e}")


def is_valid_species(animal, skip_species, breed_species):
    """
    Filter function to determine if a species should be processed

    Args:
        animal (str): Species name
        skip_species (list): Species to skip
        breed_species (list): Breed-specific species to skip

    Returns:
        bool: Whether the species should be processed
    """

    filters = [
        lambda x: not any(skip == x for skip in skip_species),
        lambda x: not any(x.startswith(breed) for breed in breed_species),
        lambda x: x == "cricetulus_griseus_crigri" if x.startswith("cricetulus_griseus_") else True,
    ]

    return all(f(animal) for f in filters)


def download_file(ftp, remote_filename, local_filename):
    if os.path.exists(local_filename):
        print(f"{local_filename} already present in {os.getcwd()}")
        return
    try:
        with open(local_filename, 'wb') as local_file:
            # Explicit RETR command with space
            ftp.retrbinary(f'RETR {remote_filename}', local_file.write)
    except ftplib.error_perm as e:
        print(f"Permission error: {e}")
        exit(1)
    except IOError as e:
        print(f"File error: {e}")
        exit(1)


def download_and_verify(ftp, filename, checksums_file):
    """
    Download file and verify its integrity

    Args:
        ftp (ftplib.FTP): FTP connection
        filename (str): Filename to download
        local_path (str): Local file path
        checksums_file (str): Checksums reference file
    """
    logging.info(f"\tDownloading: {filename}")
    download_file(ftp, filename, filename)

    # Verify checksum
    try:
        with open(checksums_file, 'r') as checksum_file:
            expected_checksum = next(
                line.split()[0]
                for line in checksum_file
                if filename in line
            )

        actual_checksum = calculate_checksum(filename)

        if actual_checksum != expected_checksum:
            print(f"Checksum mismatch for {filename}in {os.getcwd()}")
            print(f"expecte: {expected_checksum}  evaluated: {actual_checksum}")
            exit(1)

        logging.info(f"\tSuccessfully downloaded and verified: {filename}")

    except (StopIteration, ValueError) as e:
        logging.error(f"Verification failed for {filename} in {os.getcwd()}:\n{e}")
        exit(1)



def species_download(ftp, local_repository, animal, foreign_topdir):

    for dir_type in ["pep"]:  # , "dna"]:
        local_dir = os.path.join(local_repository, animal, dir_type)
        os.makedirs(local_dir, exist_ok=True)
        os.chdir(local_dir)

        foreign_dir = f"{foreign_topdir}/{animal}/{dir_type}"
        ftp.cwd(foreign_dir)

        # Download checksums
        checksum_fnm = "CHECKSUMS"
        download_file(ftp, checksum_fnm, checksum_fnm)

        # Filter and download files
        downloadable_files = list(filter(
            lambda x: x.endswith(('.dna.toplevel.gz', '.pep.all.fa.gz')), ftp.nlst()
        ))

        for item in downloadable_files:
            download_and_verify(ftp, item, checksum_fnm)


def main():
    # Configuration
    release_num = 113
    local_repository = Config.fasta_repo
    ftp_address = "ftp.ensembl.org"

    # Skip species lists
    skip_species  = Config.skip_species
    breed_species = Config.breed_species
    # Setup logging
    logging.basicConfig(
        filename='ensembl_download.log',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s: %(message)s'
    )

    # Connect to FTP
    ftp = ftp_connect(ftp_address, user='anonymous', passwd='-anonymous@')
    # move to our dir
    foreign_topdir = f"/pub/release-{release_num}/fasta"
    ftp.cwd(foreign_topdir)

    # Filter valid species using filter() and is_valid_species function
    valid_species = list(filter(
        functools.partial(is_valid_species,
                          skip_species=skip_species,
                          breed_species=breed_species),
        ftp.nlst()
    ))

    for animal in valid_species:
        print(f"\nProcessing species: {animal}")
        species_download(ftp, local_repository, animal, foreign_topdir)
        print(f"{animal} done.")

    try: ftp.quit()
    except: ftp.close()


if __name__ == "__main__":
    main()
