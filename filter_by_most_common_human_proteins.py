#!/usr/bin/env python3
import argparse
import os
import csv

# Protein list defined as a dictionary mapping the protein accession to its protein name.
# (Feel free to change this data structure if needed.)
"""
Kardell, Oliver, Thomas Gronauer, Christine von Toerne, Juliane Merl-Pham, Ann-Christine König, Teresa K. Barth, Julia Mergner, et al. “Multicenter Longitudinal Quality Assessment of MS-Based Proteomics in Plasma and Serum.” Journal of Proteome Research, February 7, 2025. https://doi.org/10.1021/acs.jproteome.4c00644.
"""
PROTEIN_LIST = {
    "O95445": "APOM",
    "P00450": "CP",
    "P00734": "F2",
    "P00738": "C1R",
    "P00747": "PLG",
    "P00748": "F12",
    "P00751": "CFB",
    "P01008": "SERPINC1",
    "P01009": "SERPINA1",
    "P01011": "SERPINA3",
    "P01019": "AGT",
    "P01023": "A2M",
    "P01024": "C3",
    "P01031": "C5",
    "P01042": "KNG1",
    "P01834": "IGKC",
    "P01857": "IGHG1",
    "P01859": "IGHG2",
    "P01860": "IGHG3",
    "P01861": "IGHG4",
    "P01871": "IGHA1",
    "P01876": "IGHM",
    "P02647": "APOA1",
    "P02649": "APOE",
    "P02652": "APOA2",
    "P02654": "APOC1",
    "P02655": "APOC2",
    "P02656": "APOC3",
    "P02743": "APCS",
    "P02749": "APOH",
    "P02750": "LRG1",
    "P02751": "IFN1",
    "P02753": "RBP4",
    "P02760": "AMBP",
    "P02763": "ORM1",
    "P02765": "AHSG",
    "P02774": "GC",
    "P02787": "ITF",
    "P02790": "HPX",
    "P03952": "KLKB1",
    "P04003": "C4BPA",
    "P04114": "APOB",
    "P04196": "HRG",
    "P04217": "A1BG",
    "P05090": "APOD",
    "P05155": "SERPING1",
    "P05158": "CFI",
    "P05546": "SERPIND1",
    "P08681": "C2",
    "P06727": "APOA4",
    "P07225": "PROS1",
    "P07357": "C8A",
    "P08185": "SERPINA6",
    "P08603": "CFH",
    "P08697": "SERPINF2",
    "P09871": "C1S",
    "P10643": "C7",
    "P10909": "CLU",
    "P13671": "C8",
    "P19852": "ORM2",
    "P19823": "ITIH2",
    "P19827": "ITIH1",
    "P20742": "PZP",
    "P22792": "CPN2",
    "P25311": "AZGP1",
    "P35822": "SERPINA4",
    "P35542": "SAA4",
    "P35858": "GFALS",
    "P43652": "AFM",
    "Q14624": "ITIH4"
}

def filter_csv_by_protein(input_csv, output_csv, protein_dict):
    """
    Reads an input CSV file and writes out only the rows where the protein accession
    (extracted from the 'ProteinID' column) is present in the protein_dict.
    
    Args:
        input_csv (str): Path to the input CSV file.
        output_csv (str): Path where the filtered CSV file will be saved.
        protein_dict (dict): Dictionary containing protein accessions to filter by.
    """
    filtered_rows = []

    # Open the input CSV for reading
    with open(input_csv, 'r', newline='') as infile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames

        for row in reader:
            # Check for either 'ProteinID' or 'Protein ID' as the column header.
            protein_id = row.get('ProteinID') or row.get('Protein ID')
            if protein_id is None:
                continue  # Skip rows that do not have a protein identifier

            # If the protein_id contains a pipe character, assume the format is like:
            # sp|ACCESSION|ENTRY and extract the accession (second field)
            if '|' in protein_id:
                parts = protein_id.split('|')
                accession = parts[1] if len(parts) > 1 else protein_id
            else:
                accession = protein_id.strip()

            # If the accession is in our protein dictionary, include this row.
            if accession in protein_dict:
                filtered_rows.append(row)

    # Write the filtered rows to the output CSV file
    with open(output_csv, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(filtered_rows)

    print(f"Filtered CSV saved to: {output_csv}")

def main():
    # Set up the argument parser with an input flag and an output directory flag.
    parser = argparse.ArgumentParser(description='Filter a CSV file by a list of proteins.')
    parser.add_argument('-i', '--input', required=True,
                        help='Path to the input CSV file.')
    parser.add_argument('-o', '--output', default='digested_glycopeptide_library',
                        help='Output directory (default: output)')
    args = parser.parse_args()

    input_csv = args.input

    # Derive the output file name by appending '_filtered' to the input filename (before extension)
    input_filename = os.path.basename(input_csv)
    name, ext = os.path.splitext(input_filename)
    output_filename = f"{name}_filtered.csv"
    output_dir = args.output

    # Create the output directory if it does not exist.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_csv = os.path.join(output_dir, output_filename)

    # Call the filtering function
    filter_csv_by_protein(input_csv, output_csv, PROTEIN_LIST)

if __name__ == '__main__':
    main()