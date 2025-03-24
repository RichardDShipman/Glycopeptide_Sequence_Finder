import os
import pandas as pd

def merge_csv_files(output_file="0_digested_glycopeptide_library.csv"):
    """Merges all CSV files in the specified directory into a single CSV file.
    
    Args:
        output_file (str, optional): Name of the output merged CSV file. Defaults to '0_digested_glycopeptide_library.csv'.
    """
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Define the relative directory where CSV files are stored
    input_directory = os.path.join(script_dir, "digested_glycopeptide_library")

    if not os.path.exists(input_directory):
        print(f"Directory '{input_directory}' does not exist.")
        return

    all_files = [f for f in os.listdir(input_directory) if f.endswith(".csv")]

    if not all_files:
        print("No CSV files found in the directory.")
        return

    merged_df = pd.DataFrame()

    for file in all_files:
        file_path = os.path.join(input_directory, file)
        df = pd.read_csv(file_path)
        df["proteome_filename_protease"] = file  # Add filename column
        merged_df = pd.concat([merged_df, df], ignore_index=True)

    output_path = os.path.join(input_directory, output_file)
    merged_df.to_csv(output_path, index=False)
    print(f"Merged CSV file saved as: {output_path}")

if __name__ == "__main__":
    merge_csv_files()