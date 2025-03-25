import argparse
import pandas as pd
import os

def load_data(glycopeptide_file, glycan_hf_file):
    """
    Load the glycopeptide data and glycan hydrophobicity data from the given files.
    """
    glycopeptide_data = pd.read_csv(glycopeptide_file, low_memory=False)
    glycan_hf_data = pd.read_csv(glycan_hf_file)
    return glycopeptide_data, glycan_hf_data

def build_glycan_hf_dict(glycan_hf_data):
    """
    Create a dictionary from glycan hydrophobicity data for fast lookup.
    """
    return dict(zip(glycan_hf_data['Glycans'], 
                    zip(glycan_hf_data['Adjusted_HF'], glycan_hf_data['Weighted_Adjusted_HF'])))

def calculate_glycopeptide_hydrophobicity(peptide_hf, glycan_composition, glycan_hf_dict, weight=1):
    """
    Calculate the combined hydrophobicity score for a glycopeptide.
    """
    if glycan_composition in glycan_hf_dict:
        glycan_hf, glycan_weighted_hf = glycan_hf_dict[glycan_composition]
        glycopeptide_hf = peptide_hf + weight * glycan_hf
        glycopeptide_weighted_hf = peptide_hf + weight * glycan_weighted_hf
        return glycopeptide_hf, glycopeptide_weighted_hf
    else:
        # Return None if glycan composition is not found
        return None, None

def calculate_scores(glycopeptide_data, glycan_hf_dict):
    """
    Calculate glycopeptide hydrophobicity scores for all glycopeptides.
    """
    glycopeptide_hf_values = []
    glycopeptide_weighted_hf_values = []

    for _, row in glycopeptide_data.iterrows():
        peptide_hf = row['Hydrophobicity']
        glycan_composition = row['composition']
        
        glycopeptide_hf, glycopeptide_weighted_hf = calculate_glycopeptide_hydrophobicity(
            peptide_hf, glycan_composition, glycan_hf_dict)
        
        glycopeptide_hf_values.append(glycopeptide_hf)
        glycopeptide_weighted_hf_values.append(glycopeptide_weighted_hf)

    return glycopeptide_hf_values, glycopeptide_weighted_hf_values

def save_results(glycopeptide_data, glycopeptide_hf_values, glycopeptide_weighted_hf_values, output_file):
    """
    Save the glycopeptide data with calculated scores to a CSV file.
    """
    glycopeptide_data['GlycopeptideHydrophobicity'] = glycopeptide_hf_values
    glycopeptide_data['GlycopeptideHydrophobicityWeighted'] = glycopeptide_weighted_hf_values
    glycopeptide_data.to_csv(output_file, index=False)
    print(f"Output saved to {output_file}")

def main():
    """
    Main function to handle the workflow of loading data, calculating scores, and saving results.
    """
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Calculate glycopeptide hydrophobicity scores.")
    parser.add_argument('-i', '--glycopeptide_file', type=str, required=True, help="Input file with glycopeptide data")
    parser.add_argument('-gh', '--glycan_hf_file', type=str, default=os.path.join(os.path.dirname(__file__), 'glycan_mass_library/human_serum_glycopeptide_hydrophobicity_data_glycan_hydrophobicity_index.csv'), help="Input file with glycan hydrophobicity data")
    parser.add_argument('-o', '--output_file', type=str, help="Output file for the results")

    args = parser.parse_args()
    
    # Load input data
    glycopeptide_data, glycan_hf_data = load_data(args.glycopeptide_file, args.glycan_hf_file)
    
    # Build glycan hydrophobicity dictionary for quick lookup
    glycan_hf_dict = build_glycan_hf_dict(glycan_hf_data)
    
    # Calculate hydrophobicity scores
    glycopeptide_hf_values, glycopeptide_weighted_hf_values = calculate_scores(glycopeptide_data, glycan_hf_dict)
    
    # Determine output file name and ensure it is saved in the HF_glycopeptide_library directory
    output_dir = os.path.join(os.path.dirname(__file__), 'HF_glycopeptide_library')
    os.makedirs(output_dir, exist_ok=True)
    output_file_name = args.output_file if args.output_file else f"{os.path.splitext(os.path.basename(args.glycopeptide_file))[0]}_HF.csv"
    output_file = os.path.join(output_dir, output_file_name)
    
    # Save the results to a CSV file
    save_results(glycopeptide_data, glycopeptide_hf_values, glycopeptide_weighted_hf_values, output_file)

if __name__ == "__main__":
    main()