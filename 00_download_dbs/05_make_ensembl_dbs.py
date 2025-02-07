#!/usr/bin/env python3

import os
import sys
import gzip
import subprocess
from time import time

from dotenv import load_dotenv
from config import Config
from el_utils.mysql import (mysql_server_connect, mysql_server_conn_close,
                            error_intolerant_search, switch_to_db, search_db, check_table_exists, count_table_rows)
from el_utils.processes import run_subprocess
from el_utils.utils import count_lines_in_compressed_file, is_gz_empty

# Load environment variables
load_dotenv()


def decompress_and_cleanup_sql(sql_gz_file):
    # Decompress SQL file
    with gzip.open(sql_gz_file, 'rb') as inf:
        sql_content = inf.read().decode()

    # Replace problematic datetime and schema references
    sql_content = sql_content.replace('0000-00-00', '1000-01-01')
    sql_content = sql_content.replace(
        'INFORMATION_SCHEMA.SESSION_VARIABLES',
        'performance_schema.session_variables'
    )

    if not sql_content:
        raise Exception(f"something went wrong when modifying {sql_gz_file}")

    return sql_content


def big_file_load(cursor, txt_file, db_name, table_name):
    came_from = os.getcwd()
    print(f"loading {txt_file} chunkwise")
    tmp_dir_name = txt_file.replace(".txt", "_chunks")
    os.makedirs(tmp_dir_name, exist_ok=True)
    os.chdir(tmp_dir_name)
    if len(os.listdir()) > 0:
        print(f"there is something in {tmp_dir_name}, and I am proceeding on the assumption that it is chunks")
    else:
        os.symlink(f"../{txt_file}", f"./{txt_file}")
        cmd = f"split -C 1G {txt_file}"
        time0 = time()
        run_subprocess(cmd)
        print(f"split {txt_file} in 1G chunks in {time()-time0:.0f} secs.")
        os.unlink(f"./{txt_file}")
    dirs = sorted(list(os.listdir()))
    start_idx = 0
    start_idx = dirs.index('xab') + 1

    for chunk in dirs[start_idx:]:
        time0 = time()
        print(chunk)
        fullpath = f"{os.getcwd()}/{chunk}"
        load_qry = f"LOAD DATA LOCAL INFILE '{fullpath}' INTO TABLE {db_name}.{table_name}"
        search_db(cursor, load_qry)  # we'll fail dn duplicates
        cursor.connection.commit()
        print(f"loaded chunk {chunk} in {time()-time0:.0f} secs.")

    os.chdir(came_from)


def load_data(cursor, db_name, txt_gz_file, dry_run):

    # if is_gz_empty(txt_gz_file):
    #     print(f"{txt_gz_file} seems to be empty.")
    #     return

    txt_file   = txt_gz_file.replace('.gz', '')
    fullpath   = f"{os.getcwd()}/{txt_file}"
    table_name = txt_file.replace('.txt', '')
    if not check_table_exists(cursor, db_name, table_name):
        print(f"table {table_name} not found in {db_name}")
        # exit()
        return  # it is not our job here to make the table

    load_qry = f"LOAD DATA LOCAL INFILE '{fullpath}' INTO TABLE {db_name}.{table_name}"
    if dry_run:
        print(f"decompressing and loading {fullpath}")
        print(load_qry)
        return

    # check first if the table is maybe loaded already
    if (number_of_rows_in_db := count_table_rows(cursor, db_name, table_name)) > 0:
        print(f"I am counting lines in {txt_gz_file} ...")
        number_of_lines_in_file = count_lines_in_compressed_file(txt_gz_file)
        if number_of_rows_in_db == number_of_lines_in_file:
            print(f"Table {table_name} ok. The number of rows is {number_of_rows_in_db}.")
            return
        else:
            # delete data from the table, because the alternative is to check row-by_row
            print(f"Table {table_name} was not loaded properly.")
            print(f"Number of rows in the table {number_of_rows_in_db}. Lines in the file {number_of_lines_in_file}.")
            print(f"Deleting the rows from the table, and attempting the re-load.")
            error_intolerant_search(cursor, f"delete from {db_name}.{table_name}")
    else:
        print(f"table {table_name} exists in {db_name}, but is empty")

    # Decompress text file
    cmd = f"gunzip -c {txt_gz_file}"
    run_subprocess(cmd, stdoutfnm=txt_file, noexit=True)

    # What is the size of the decompressed file?
    if not os.path.exists(txt_file):
        print(f"error gunzipping {txt_gz_file} in {os.getcwd()} - gunzipped file not produced:\n{cmd}")
        exit()

    decompressed_size = os.path.getsize(txt_file)
    if decompressed_size == 0:
        print(f"error gunzipping {txt_gz_file} in {os.getcwd()} - gunzipped file is empty:\n{cmd}")
        exit()

    size_in_g = decompressed_size / 1024**3
    if size_in_g < 1:  # less than a gigabyte
        print(f"size of {txt_file} is {size_in_g}G; loading directly")
        search_db(cursor, load_qry)
    else:
        print(f"size of {txt_file} is {size_in_g}G; will receive special treatment")
        big_file_load(cursor, txt_file, db_name, table_name)
    # Remove temporary text file
    if os.path.exists(txt_file): os.remove(txt_file)


def get_sql_file(db):
    sql_files = [fnm for fnm in os.listdir() if fnm.startswith(db) and fnm.endswith("sql.gz")]
    if len(sql_files) == 0:
        print(f"no sql file found in {os.getcwd()}")
        exit(1)
    elif len(sql_files) > 1:
        print(f"multiple sql files found in {os.getcwd()}: {sql_files}")
        exit(1)
    return sql_files[0]


def make_tables(cursor, sql_gz_file, dry_run):
    sql_content = ""
    if dry_run:
        print("decompress and cleanup the sql file")
    else:
        sql_content = decompress_and_cleanup_sql(sql_gz_file)

    # Import SQL file
    if dry_run:
        print(f"execute the contents of the sql file - create tables")
    else:
        error_intolerant_search(cursor, sql_content)


def main():

    dry_run = False
    from_scratch = False  # if this is True, we will drop the exiting tables and make new ones

    path_to_db  = Config.mysql_repo
    # Validate database path exists
    if not os.path.exists(path_to_db):
        raise FileNotFoundError(f"{path_to_db} not found")
    os.chdir(path_to_db)
    # Get list of databases
    dbs = sorted([d for d in os.listdir('.') if os.path.isdir(d)])

    # MySQL connection parameters from .env
    cursor = mysql_server_connect(user=os.getenv('MYSQL_USER'),
                                  passwd=os.getenv('MYSQL_PASSWORD'),
                                  host=os.getenv('MYSQL_HOST', 'localhost'),
                                  port=int(os.getenv('MYSQL_PORT', 3306)))

    error_intolerant_search(cursor, "set GLOBAL local_infile = 'ON'")  # allow loading from non-privileged dir
    # the following are supposed to make the loading faster
    error_intolerant_search(cursor, "set autocommit = 0")  # allow loading from non-privileged dir
    print("autocommit", error_intolerant_search(cursor, "select @@autocommit")[0][0])

    error_intolerant_search(cursor, "set unique_checks=0")
    print("unique_checks", error_intolerant_search(cursor, "select @@unique_checks")[0][0])

    error_intolerant_search(cursor, "set foreign_key_checks = 0")
    print("foreign_key_checks", error_intolerant_search(cursor, "select @@foreign_key_checks")[0][0])

    for db in dbs:
        if not db.startswith("ensembl"): continue
        print()
        print(f"************************")
        print(db)

        os.chdir(os.path.join(path_to_db, db))
        # find the SQL file exists
        sql_gz_file = get_sql_file(db)
        db_name = sql_gz_file.replace(".sql.gz", "")
        # Create database
        qry = f"CREATE DATABASE IF NOT EXISTS `{db_name}`"
        if dry_run:
            print(qry)
        else:
            error_intolerant_search(cursor, qry)
            switch_to_db(cursor, db_name)

        if from_scratch:
            make_tables(cursor, sql_gz_file, dry_run)

        # Process and import text files
        txt_files = [f for f in os.listdir('.') if f.endswith('.txt.gz')]

        for txt_gz_file in txt_files:
            # TODO see if I can drop the cigar line and pctg as floats to make homology member smaller
            # TODO switch to some kind of ORM  and move to postgres
            if "homology" in txt_gz_file: continue  # big, the loading is a clunker
            print()
            print(f"loading {txt_gz_file}")
            load_data(cursor, db_name, txt_gz_file, dry_run)
            print(f"loading {txt_gz_file} done")

        print(f"{db} done")

    mysql_server_conn_close(cursor)


if __name__ == "__main__":
    main()
