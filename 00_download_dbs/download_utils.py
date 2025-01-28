import ftplib
import logging
import os
import subprocess


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
        lambda x: not x.startswith("ensembl_"),
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
