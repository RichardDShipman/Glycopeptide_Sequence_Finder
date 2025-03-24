import os
import sqlite3
import argparse
import pandas as pd

# Function to merge CSV files into an SQLite database
def merge_csv_to_sqlite(directory, output_db):
    # Create SQLite connection and cursor
    conn = sqlite3.connect(output_db)
    cursor = conn.cursor()

    # Create a master table to track filenames and associated metadata
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS master_table (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        filename TEXT NOT NULL,
        row_count INTEGER NOT NULL,
        column_count INTEGER NOT NULL
    )
    """)

    # Create a table to store all the data
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS glycopeptide_data (
        ProteinID TEXT,
        Site REAL,
        GlyToucan_AC TEXT,
        Composition TEXT,
        ShorthandGlycan TEXT,
        Peptide TEXT,
        Start REAL,
        End REAL,
        Length REAL,
        Sequon TEXT,
        GlycopeptideMass REAL,
        PeptideMass REAL,
        GlycanMass REAL,
        Hydrophobicity REAL,
        pI REAL,
        z2 REAL,
        Charge REAL,
        IonSeries TEXT
    )
    """)

    # Process each CSV file in the directory
    for file in os.listdir(directory):
        if file.endswith('.csv'):
            file_path = os.path.join(directory, file)
            print(f"Processing {file_path}")

            # Read the CSV file into a DataFrame
            df = pd.read_csv(file_path, encoding='ISO-8859-1')

            # Insert metadata into the master table
            cursor.execute("""
            INSERT INTO master_table (filename, row_count, column_count)
            VALUES (?, ?, ?)
            """, (file, len(df), len(df.columns)))

            # Insert data into the glycopeptide_data table
            df.to_sql('glycopeptide_data', conn, if_exists='append', index=False)

    # Commit and close the connection
    conn.commit()
    conn.close()
    print(f"All files have been processed and merged into {output_db}")

# Main function to parse command-line arguments and execute the script
def main():
    # Set up argument parser for command-line inputs
    parser = argparse.ArgumentParser(description="Merge CSV files into an SQLite database.")
    parser.add_argument("-i", "--input_dir", type=str, default="digested_glycopeptide_library", help="Directory containing CSV files.")
    parser.add_argument("-o", "--output_db", type=str, default=None, help="Output SQLite database file. Defaults to '0_digested_glycopeptide_library_sqlite.db'.")

    # Parse the arguments
    args = parser.parse_args()

    # Set output database name if not provided
    if not args.output_db:
        args.output_db = os.path.join(args.input_dir, "0_digested_glycopeptide_library_sqlite.db")

    # Ensure the input directory exists
    if not os.path.isdir(args.input_dir):
        print(f"Error: The directory '{args.input_dir}' does not exist.")
        return

    # Call the merge function
    merge_csv_to_sqlite(args.input_dir, args.output_db)

if __name__ == "__main__":
    main()
