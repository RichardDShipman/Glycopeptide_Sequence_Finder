import os
import sqlite3
import argparse
import pandas as pd
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def merge_csv_to_sqlite(glycopeptide_dir, glycosite_dir, output_db):
    # Create SQLite connection and cursor
    conn = sqlite3.connect(output_db)
    cursor = conn.cursor()

    # Define the schema for the glycopeptide_data table
    glycopeptide_columns = [
        "ProteinID", "Site", "GlyToucan_AC", "Composition", "ShorthandGlycan",
        "Peptide", "Start", "End", "Length", "Sequon", "GlycopeptideMass",
        "PeptideMass", "GlycanMass", "Hydrophobicity", "pI", "z2", "Charge", "IonSeries"
    ]

    # Define the schema for the glycosite_data table
    glycosite_columns = [
        "ProteinID", "Site", "Sequon", "Species", "TaxonID", "GeneName", "seq21"
    ]

    # Drop existing tables if they exist
    cursor.execute("DROP TABLE IF EXISTS master_table")
    cursor.execute("DROP TABLE IF EXISTS glycopeptide_data")
    cursor.execute("DROP TABLE IF EXISTS glycosite_data")

    # Create a master table to track filenames and associated metadata
    cursor.execute("""
    CREATE TABLE master_table (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        filename TEXT NOT NULL,
        row_count INTEGER NOT NULL,
        column_count INTEGER NOT NULL,
        status TEXT NOT NULL,
        error_message TEXT,
        table_name TEXT NOT NULL
    )
    """)

    # Create tables to store the data
    cursor.execute(f"""
    CREATE TABLE glycopeptide_data (
        {', '.join([f'{col} REAL' if col != 'ProteinID' and col != 'GlyToucan_AC' and col != 'Composition' and col != 'ShorthandGlycan' and col != 'Peptide' and col != 'Sequon' and col != 'IonSeries' else f'{col} TEXT' for col in glycopeptide_columns])}
    )
    """)

    cursor.execute(f"""
    CREATE TABLE glycosite_data (
        {', '.join([f'{col} TEXT' for col in glycosite_columns])}
    )
    """)

    # Process glycopeptide CSV files
    process_directory(glycopeptide_dir, cursor, conn, glycopeptide_columns, 'glycopeptide_data')

    # Process glycosite CSV files
    process_directory(glycosite_dir, cursor, conn, glycosite_columns, 'glycosite_data')

    # Commit and close the connection
    conn.commit()
    conn.close()
    logging.info(f"All files have been processed and merged into {output_db}")

def process_directory(directory, cursor, conn, columns, table_name):
    """Process all CSV files in a directory and insert them into the specified table."""
    if not os.path.isdir(directory):
        logging.warning(f"Directory does not exist: {directory}")
        return

    for file in os.listdir(directory):
        if file.endswith('.csv'):
            file_path = os.path.join(directory, file)
            logging.info(f"Processing {file_path}")

            try:
                # Check if file is empty
                if os.stat(file_path).st_size == 0:
                    logging.warning(f"Empty file encountered: {file_path}")
                    cursor.execute("""
                    INSERT INTO master_table (filename, row_count, column_count, status, error_message, table_name)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """, (file, 0, 0, "EMPTY", "File is empty", table_name))
                    continue

                # Read the CSV file into a DataFrame
                df = pd.read_csv(file_path, encoding='ISO-8859-1')

                # Check if DataFrame is empty
                if df.empty:
                    logging.warning(f"Empty DataFrame from file: {file_path}")
                    cursor.execute("""
                    INSERT INTO master_table (filename, row_count, column_count, status, error_message, table_name)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """, (file, 0, 0, "EMPTY_DF", "DataFrame is empty", table_name))
                    continue

                # Align DataFrame columns with the table schema
                df = df.reindex(columns=columns, fill_value=None)

                # Insert metadata into the master table
                cursor.execute("""
                INSERT INTO master_table (filename, row_count, column_count, status, error_message, table_name)
                VALUES (?, ?, ?, ?, ?, ?)
                """, (file, len(df), len(df.columns), "SUCCESS", None, table_name))

                # Insert data into the appropriate table
                df.to_sql(table_name, conn, if_exists='append', index=False)
                logging.info(f"Successfully processed {file_path}: {len(df)} rows")

            except pd.errors.EmptyDataError:
                logging.warning(f"EmptyDataError for file: {file_path}")
                cursor.execute("""
                INSERT INTO master_table (filename, row_count, column_count, status, error_message, table_name)
                VALUES (?, ?, ?, ?, ?, ?)
                """, (file, 0, 0, "ERROR", "EmptyDataError: No columns to parse", table_name))
            except Exception as e:
                logging.error(f"Error processing {file_path}: {str(e)}")
                cursor.execute("""
                INSERT INTO master_table (filename, row_count, column_count, status, error_message, table_name)
                VALUES (?, ?, ?, ?, ?, ?)
                """, (file, 0, 0, "ERROR", str(e), table_name))

def main():
    # Set up argument parser for command-line inputs
    parser = argparse.ArgumentParser(description="Merge CSV files into an SQLite database.")
    parser.add_argument("-i", "--input_dir", type=str, default="./data/digested_glycopeptide_library", help="Directory containing glycopeptide CSV files.")
    parser.add_argument("-g", "--glycosite_dir", type=str, default="./data/glycosite_library", help="Directory containing glycosite CSV files.")
    parser.add_argument("-o", "--output_db", type=str, default=None, help="Output SQLite database file. Defaults to '0_digested_glycopeptide_library_sqlite.db'.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Enable verbose logging.")

    # Parse the arguments
    args = parser.parse_args()

    # Set logging level based on verbose flag
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Set output database name if not provided
    if not args.output_db:
        args.output_db = os.path.join(args.input_dir, "0_digested_glycopeptide_library_sqlite.db")

    # Call the merge function
    merge_csv_to_sqlite(args.input_dir, args.glycosite_dir, args.output_db)

if __name__ == "__main__":
    main()
