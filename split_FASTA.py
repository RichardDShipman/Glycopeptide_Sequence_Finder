import argparse
import os
import math

def split_fasta(input_file, splits):
    """
    Split a FASTA file into smaller parts.

    Parameters:
    - input_file: str, the path to the input protein FASTA file.
    - splits: int, the number of parts to split the FASTA file into.

    Returns:
    - None, saves the resulting splits as separate FASTA files.
    """
    # Read the entire file into memory
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Ensure the file contains sequence data
    if len(lines) == 0 or not any(line.startswith('>') for line in lines):
        raise ValueError("Invalid FASTA file format.")

    # Determine the number of records (sequences) in the FASTA file
    records = []
    current_record = []
    for line in lines:
        if line.startswith('>'):
            if current_record:
                records.append(current_record)
            current_record = [line]
        else:
            current_record.append(line)
    if current_record:
        records.append(current_record)

    # Split the records into approximately equal parts
    total_records = len(records)
    records_per_split = math.ceil(total_records / splits)

    # Write the splits into separate files
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    for i in range(splits):
        start = i * records_per_split
        end = min((i + 1) * records_per_split, total_records)
        split_file = f"{base_name}_part_{i + 1}.fasta"

        with open(split_file, 'w') as out_f:
            for record in records[start:end]:
                for line in record:
                    out_f.write(line)

        print(f"Created {split_file} with {end - start} records.")

def main():
    """
    Main function to parse arguments and execute the splitting function.
    """
    parser = argparse.ArgumentParser(description="Split a protein FASTA file into smaller files.")
    parser.add_argument('-i', '--input', type=str, required=True, help="Input protein FASTA file.")
    parser.add_argument('-s', '--splits', type=int, default=2, help="Number of splits to create (default: 2).")
    
    args = parser.parse_args()

    # Call the split function with provided arguments
    try:
        split_fasta(args.input, args.splits)
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
