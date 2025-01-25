#!/usr/bin/env python3
# this is untested;  move to mysqldn and implmentend stufff
import os
import sys
import gzip
import subprocess
import mysql.connector
from dotenv import load_dotenv

# Load environment variables
load_dotenv()


def main():
    release_num = 110
    path_to_db = f"/storage/databases/ensembl-{release_num}/mysql"

    # Validate database path exists
    if not os.path.exists(path_to_db):
        raise FileNotFoundError(f"{path_to_db} not found")

    os.chdir(path_to_db)

    dry_run = False

    # Get list of databases
    dbs = [d for d in os.listdir('.') if os.path.isdir(d)]

    # MySQL connection parameters from .env
    db_config = {
        'user': os.getenv('MYSQL_USER'),
        'password': os.getenv('MYSQL_PASSWORD'),
        'host': os.getenv('MYSQL_HOST', 'localhost')
    }

    for db in dbs:
        print(f"************************")
        print(db)

        os.chdir(os.path.join(path_to_db, db))

        # Create database
        if not dry_run:
            try:
                with mysql.connector.connect(**db_config) as connection:
                    with connection.cursor() as cursor:
                        cursor.execute(f"CREATE DATABASE IF NOT EXISTS `{db}`")
            except mysql.connector.Error as err:
                print(f"Error creating database {db}: {err}")
                continue

        # Decompress and process SQL file
        sql_gz_file = f"{db}.sql.gz"
        if not os.path.exists(sql_gz_file):
            raise FileNotFoundError(f"{sql_gz_file} not found")

        # Decompress SQL file
        with gzip.open(sql_gz_file, 'rb') as f_in:
            with open(f"{db}.sql", 'wb') as f_out:
                f_out.write(f_in.read())

        # Modify SQL file contents
        with open(f"{db}.sql", 'r') as f:
            sql_content = f.read()

        # Replace problematic datetime and schema references
        sql_content = sql_content.replace('0000-00-00', '1000-01-01')
        sql_content = sql_content.replace(
            'INFORMATION_SCHEMA.SESSION_VARIABLES',
            'performance_schema.session_variables'
        )

        # Import SQL file
        if not dry_run:
            try:
                with mysql.connector.connect(**db_config, database=db) as connection:
                    with connection.cursor() as cursor:
                        cursor.execute(sql_content)
            except mysql.connector.Error as err:
                print(f"Error importing SQL for {db}: {err}")

        # Remove temporary SQL file
        os.remove(f"{db}.sql")

        # Process and import text files
        txt_files = [f for f in os.listdir('.') if f.endswith('.txt.gz')]

        for txt_gz_file in txt_files:
            # Decompress text file
            with gzip.open(txt_gz_file, 'rb') as f_in:
                txt_file = txt_gz_file.replace('.gz', '')
                with open(txt_file, 'wb') as f_out:
                    f_out.write(f_in.read())

            # Import text file using mysqlimport equivalent
            if not dry_run:
                try:
                    with mysql.connector.connect(**db_config, database=db) as connection:
                        with connection.cursor() as cursor:
                            with open(txt_file, 'r') as f:
                                for line in f:
                                    # Basic import - you might need to adjust based on exact requirements
                                    cursor.execute(
                                        f"LOAD DATA LOCAL INFILE '{txt_file}' INTO TABLE {txt_file.replace('.txt', '')}")
                except mysql.connector.Error as err:
                    print(f"Error importing text file {txt_file}: {err}")

            # Remove temporary text file
            os.remove(txt_file)

        print()


if __name__ == "__main__":
    main()
