import numpy as np
import pandas as pd
import re
import argparse
import os

# Define Encoding Definitions

# Define Peptide (amino acid dictionary) (20 standard amino acids)
AA_DICT = {aa: i for i, aa in enumerate("ACDEFGHIKLMNPQRSTVWY")} # 'B' (Asx), 'Z' (Glx), 'U' (Sec), 'O' (Pyl) removed for simplicity
MAX_PEPTIDE_LENGTH = 50 # defines the maximum length of the peptide sequence

# Define Glycan (monosaccharide dictionary)
GLYCAN_MONOSACCHARIDES = ['N', 'H', 'F', 'A'] # 'G' (NeuGc), 'X' (Pent/Xylo), 'Ph' (Phospho), 'Su' (Sulfo), 'aH' (Hexuronic acid) removed for simplicity
MAX_GLYCAN_LENGTH = 20 # defines the maximum length of the glycan sequence

# Charge states (Max charge = 5)
MAX_CHARGE = 5 # defines max charge state in the encoding

def encode_peptide(peptide):
    """Encode peptide composition as a 50x20 matrix."""
    
    # Initialize the encoding matrix
    encoding = np.zeros((MAX_PEPTIDE_LENGTH, len(AA_DICT)))
    
    # Encode each amino acid in the peptide sequence
    for i, aa in enumerate(peptide[:MAX_PEPTIDE_LENGTH]):
        if aa in AA_DICT:
            encoding[i, AA_DICT[aa]] = 1
    return encoding

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

def encode_glycan(glycan_sequence):
    # Define the possible monosaccharides and their corresponding index
    monosaccharides = GLYCAN_MONOSACCHARIDES
    encoding = np.zeros((MAX_GLYCAN_LENGTH, len(GLYCAN_MONOSACCHARIDES)), dtype=int)  # Initialize the encoding matrix
    
    # Iterate through the sequence and encode each monosaccharide
    position = 0
    for i, monosaccharide in enumerate(glycan_sequence):
        if position >= MAX_GLYCAN_LENGTH:
            break  # Stop encoding if we exceed the maximum length
        # Set the correct column for the current monosaccharide
        if monosaccharide in monosaccharides:
            # The index of the monosaccharide in the row (N=0, H=1, F=2, A=3)
            monosaccharide_idx = monosaccharides.index(monosaccharide)
            encoding[position, monosaccharide_idx] = 1
            position += 1
    
    return encoding

def encode_charge(charge):
    """Encode charge state (Z) as a one-hot vector with 10 columns."""
    encoding = np.zeros((1, MAX_CHARGE))
    if 1 <= charge <= MAX_CHARGE:
        encoding[0, charge - 1] = 1
    return encoding

def encode_glycopeptide(peptide, glycan, charge):
    """Create a feature vector from peptide composition, glycan composition, and charge state."""
    
    # Encode peptide, glycan, and charge
    peptide_encoded = encode_peptide(peptide)
    glycan_encoded = encode_glycan(glycan)
    charge_encoded = encode_charge(charge)

    # Concatenate the encoded features into a single vector
    feature_vector = np.concatenate([peptide_encoded.flatten(), glycan_encoded.flatten(), charge_encoded.flatten()])
    
    return feature_vector

def generate_encoding_definition(output_file):
    """Generate a CSV file that defines the encoding format."""
    encoding_definitions = []

    # Peptide encoding
    for aa, idx in AA_DICT.items():
        encoding_definitions.append({'Type': 'Peptide', 'Position': f'1-{MAX_PEPTIDE_LENGTH}', 'Feature': aa, 'Index': idx, 'Size': len(AA_DICT)})

    # Glycan encoding
    for idx, mono in enumerate(GLYCAN_MONOSACCHARIDES):
        encoding_definitions.append({'Type': 'Glycan', 'Position': f'1-{MAX_GLYCAN_LENGTH}', 'Feature': mono, 'Index': idx, 'Size': len(GLYCAN_MONOSACCHARIDES)})

    # Charge encoding
    for i in range(1, MAX_CHARGE + 1):  # Charge states from 1 to MAX_CHARGE
        encoding_definitions.append({'Type': 'Charge', 'Position': '1', 'Feature': f'Charge_{i}', 'Index': i-1, 'Size': MAX_CHARGE})

    # Add encoding size row
    encoding_size = len(encoding_definitions)
    encoding_definitions.append({'Type': 'Encoding_Size', 'Position': '-', 'Feature': '-', 'Index': '-', 'Size': encoding_size})

    # Convert to DataFrame and save as CSV
    df_definitions = pd.DataFrame(encoding_definitions)
    df_definitions.to_csv(output_file, index=False)
    print(f"Encoding definition saved as '{output_file}', total encoding size: {encoding_size}")

def count_ones_zeros(one_hot_vector):
    """Count the number of 1s and 0s in a one-hot encoded vector."""
    ones_count = np.sum(one_hot_vector)  # Count total 1s
    zeros_count = len(one_hot_vector) - ones_count  # Total size - 1s gives 0s
    return ones_count, zeros_count

def main():
    parser = argparse.ArgumentParser(description="Encode glycopeptides with one-hot encoding.")
    parser.add_argument('-i', '--input', required=True, help="Input CSV file")
    parser.add_argument('-o', '--output', help="Output CSV file (default: input file name with '_onehotencodings.csv' in 'one_hot_encodings' directory)")
    parser.add_argument('-p', '--peptide', default='Peptide', help="Column name for peptide sequence (default: 'peptide')")
    parser.add_argument('-g', '--glycan', default='ShorthandGlycan', help="Column name for glycan composition (default: 'glycan')")
    parser.add_argument('-z', '--charge', default='Charge', help="Column name for charge state (default: 'charge')")
    parser.add_argument('-d', '--definition', nargs='?', const=os.path.join('one_hot_encodings', 'encoding_definition.txt'), help="Output CSV file for encoding definition (default: 'one_hot_encodings/encoding_definition.txt')")

    # Parse arguments
    args = parser.parse_args()
    
    input_file = args.input
    output_dir = 'one_hot_encodings'
    os.makedirs(output_dir, exist_ok=True)
    output_file = args.output if args.output else os.path.join(output_dir, os.path.splitext(os.path.basename(input_file))[0] + '_onehotencoded.csv')
    peptide_col = args.peptide
    glycan_col = args.glycan
    charge_col = args.charge
    
    # Read input data
    df = pd.read_csv(input_file)
    
    # Encode all glycopeptides
    encoded_data = np.array([encode_glycopeptide(row[peptide_col], glycan_composition_to_sequence(row[glycan_col]), row[charge_col]) for _, row in df.iterrows()])
    
    # Add glycan composition sequence as a new column
    df['Glycan_Composition_Sequence'] = df[glycan_col].apply(glycan_composition_to_sequence)

    # Add glycopeptide sequence + charge with padding
    df['Peptide_Padded'] = df[peptide_col].apply(lambda x: x.ljust(MAX_PEPTIDE_LENGTH, 'X'))
    df['Glycan_Padded'] = df['Glycan_Composition_Sequence'].apply(lambda x: x.ljust(MAX_GLYCAN_LENGTH, 'X'))
    df['Glycopeptide_Composition_Sequence'] = df['Peptide_Padded'] + df['Glycan_Padded'] + df[charge_col].astype(str)

    # Add the one-hot encoded data as a new column
    df['One_Hot_Encoding'] = [list(encoded_data[i]) for i in range(len(df))]

    # Add the number of 1s and 0s in the one-hot encoding
    df['Ones_Count'] = df['One_Hot_Encoding'].apply(lambda x: count_ones_zeros(x)[0])
    df['Zeros_Count'] = df['One_Hot_Encoding'].apply(lambda x: count_ones_zeros(x)[1])
    
    # Save the output
    df.to_csv(output_file, index=False)
    print(f"Encoded data saved as '{output_file}'")

    # Generate encoding definition if requested
    if args.definition:
        generate_encoding_definition(args.definition)

if __name__ == "__main__":
    main()
