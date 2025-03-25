import pandas as pd
import argparse
import os
from scipy.stats import zscore

def compute_adjusted_hf(df):
    """Compute the adjusted HF score by removing the peptide effect."""
    peptide_means = df.groupby('Peptide Backbone')['HFs'].transform('mean')
    df['Adjusted_HF'] = df['HFs'] - peptide_means
    return df

def compute_weighted_adjusted_hf(df):
    """Compute weighted adjusted HF using frequency-based adjustment and Z-score normalization."""
    df['Z_Adjusted_HF'] = zscore(df['Adjusted_HF'])  # Normalize adjusted HF values
    glycan_counts = df['Glycans'].value_counts().to_dict()
    df['Glycan_Frequency'] = df['Glycans'].map(glycan_counts)
    df['Weighted_Adjusted_HF'] = df['Z_Adjusted_HF'] * df['Glycan_Frequency']
    # calcuate a scaled version of the weighted adjusted HF between -1 and 1
    df['Scaled_HF'] = df['Weighted_Adjusted_HF'] / df['Weighted_Adjusted_HF'].abs().max()

    return df

def count_glycan_frequency(df):
    """Count the frequency of each glycan in the dataset."""
    glycan_counts = df['Glycans'].value_counts().reset_index()
    glycan_counts.columns = ['Glycans', 'Count']
    return glycan_counts

def rank_glycans_by_hydrophobicity(input_csv, output_csv=None):
    """Rank glycans based on their adjusted hydrophobicity factor."""
    df = pd.read_csv(input_csv)
    
    # Compute the adjusted HF score by removing peptide effect
    df = compute_adjusted_hf(df)
    
    # Compute the weighted adjusted HF
    df = compute_weighted_adjusted_hf(df)
    
    # Group glycans by structure and compute mean adjusted HF and weighted adjusted HF
    glycan_summary = df.groupby('Glycans').agg(
        Adjusted_HF=('Adjusted_HF', 'mean'),
        Z_Adjusted_HF=('Z_Adjusted_HF', 'mean'),
        Weighted_Adjusted_HF=('Weighted_Adjusted_HF', 'mean'),
        Scaled_HF=('Scaled_HF', 'mean')
    ).reset_index()
    
    # Count glycan occurrences
    glycan_counts = count_glycan_frequency(df)
    
    # Merge count data with glycan summary
    glycan_summary = glycan_summary.merge(glycan_counts, on='Glycans', how='left')
    
    # Sort glycans by weighted adjusted HF score (ascending)
    glycan_summary = glycan_summary.sort_values(by='Scaled_HF', ascending=True)
    
    # Determine the output file name if not provided
    if output_csv is None:
        base = os.path.splitext(input_csv)[0]
        output_csv = f"{base}_glycan_hydrophobicity_index.csv"
    
    # Save the results to a new CSV file
    glycan_summary.to_csv(output_csv, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Rank glycans by hydrophobicity.')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file path')
    parser.add_argument('-o', '--output', help='Output CSV file path')
    args = parser.parse_args()
    
    rank_glycans_by_hydrophobicity(args.input, args.output)
