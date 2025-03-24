import csv
import argparse
import os

# Function to split the 'byonic' column into 'composition' and 'mass' columns
def split_byonic(input_file, output_file):
    """Splits the 'byonic' column into 'composition' and 'mass' columns in a CSV file."""
    
    # Open the input and output CSV files
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)  # Read the input CSV file as a dictionary
        fieldnames = reader.fieldnames + ['composition', 'mass']  # Add new fieldnames for the output CSV
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)  # Create a writer object with the new fieldnames
        
        writer.writeheader()  # Write the header to the output CSV file
        for row in reader:
            byonic = row.get('byonic') or row.get('sequence_byonic')  # Get the 'byonic' or 'sequence_byonic' column value
            if '%' in byonic:
                composition, mass = byonic.split('%')  # Split the 'byonic' value into 'composition' and 'mass'
                row['composition'] = composition.strip()  # Strip any extra spaces and assign to 'composition'
                row['mass'] = mass.strip()  # Strip any extra spaces and assign to 'mass'
            else:
                row['composition'] = ''  # If no '%' is found, set 'composition' to empty
                row['mass'] = ''  # If no '%' is found, set 'mass' to empty
            writer.writerow(row)  # Write the modified row to the output CSV file

# Main function to handle command-line arguments and call the split_byonic function
def main():
    parser = argparse.ArgumentParser(description='Split byonic column into composition and mass columns.')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file')  # Input file argument
    parser.add_argument('-o', '--output', help='Output CSV file')  # Output file argument
    args = parser.parse_args()  # Parse the command-line arguments
    
    # If output file is not provided, create a default output file name
    if not args.output:
        input_filename = os.path.splitext(os.path.basename(args.input))[0]
        args.output = f"{input_filename}_glycan_library.csv"
    
    split_byonic(args.input, args.output)  # Call the split_byonic function with the provided arguments
    print(f'Output written to: {args.output}')  # Print the output file path

# Entry point of the script
if __name__ == "__main__":
    main()