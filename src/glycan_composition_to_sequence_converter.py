import argparse
import pandas as pd
import re

GLYCAN_MONOSACCHARIDES = ['N', 'H', 'F', 'A', 'S', 'G', 'X', 'Ph', 'Su', 'aH'] 

def glycan_composition_to_sequence(glycan_composition):
    # Define the mapping for the monosaccharides
    monosaccharides = GLYCAN_MONOSACCHARIDES
    
    # Create a regex pattern based on the monosaccharides list
    pattern = '|'.join(monosaccharides)
    
    # Initialize the glycan sequence dictionary
    glycan_sequence_dict = {mono: 0 for mono in monosaccharides}
    
    # Find all monosaccharide components and their counts using regex
    components = re.findall(rf'({pattern})(\d+)', glycan_composition)
    
    for component in components:
        monosaccharide = component[0]
        count = int(component[1])
        
        # Update the count for the monosaccharide in the dictionary
        glycan_sequence_dict[monosaccharide] += count
    
    # Create the ordered glycan sequence based on GLYCAN_MONOSACCHARIDES
    glycan_sequence = ''.join([mono * glycan_sequence_dict[mono] for mono in monosaccharides])
    
    return glycan_sequence

def process_glycan_data(input_file, output_file, glycan_column):
    # Load the CSV file into a pandas DataFrame
    df = pd.read_csv(input_file)
    
    # Check if the specified glycan column exists
    if glycan_column not in df.columns:
        raise ValueError(f"Column '{glycan_column}' not found in the input file.")
    
    # Apply the glycan composition to sequence function
    df['glycan_composition_sequence'] = df[glycan_column].apply(glycan_composition_to_sequence)
    
    # Save the output to a new CSV file
    df.to_csv(output_file, index=False)

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Convert glycan composition to glycan sequence.')
    
    # Add arguments for input, output, and glycan column
    parser.add_argument('-i', '--input', required=True, help='Input CSV file path')
    parser.add_argument('-o', '--output', default=None, help='Output CSV file path (default: input filename + "_composition_sequences.csv")')
    parser.add_argument('-g', '--glycan', required=True, help='Column name containing glycan composition data')
    
    # Parse arguments
    args = parser.parse_args()

    # Set the default output file if not provided
    if not args.output:
        output_file = f"{args.input.split('.')[0]}_composition_sequences.csv"
    else:
        output_file = args.output
    
    # Process the glycan data
    process_glycan_data(args.input, output_file, args.glycan)
    print(f"Converted glycan data saved to {output_file}")

if __name__ == '__main__':
    main()