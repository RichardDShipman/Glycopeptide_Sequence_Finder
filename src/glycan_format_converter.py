import pandas as pd
import argparse
import re
import os

# Mapping definitions
SHORTHAND_TO_LONGHAND = {
    'N': 'HexNAc',
    'H': 'Hex',
    'F': 'Fuc',
    'A': 'NeuAc',
    'G': 'NeuGc',
    'X': 'Pent',
    'Ph': 'Phospho',
    'Su': 'Sulpho',
    'aH': 'HexA'
}

# Reverse mapping
LONGHAND_TO_SHORTHAND = {v: k for k, v in SHORTHAND_TO_LONGHAND.items()}

def convert_shorthand_to_longhand(shorthand, order=None):
    """Converts shorthand glycan notation (e.g., H5N4A1) to longhand (e.g., HexNAc(4)Hex(5)NeuAc(1))."""
    matches = re.findall(r'([A-Za-z]+)(\d*)', shorthand)
    default_order = ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc', 'Pent', 'Phospho', 'Sulpho']
    order = order or default_order
    longhand_parts = []
    for code, count in matches:
        if code in SHORTHAND_TO_LONGHAND:
            longhand_parts.append((SHORTHAND_TO_LONGHAND[code], count or '1'))
    ordered_longhand_parts = sorted(longhand_parts, key=lambda x: order.index(x[0]) if x[0] in order else len(order))
    return ''.join(f"{name}({count})" for name, count in ordered_longhand_parts)

def convert_longhand_to_shorthand(longhand, order=None):
    """Converts longhand glycan notation (e.g., HexNAc(4)Hex(5)NeuAc(1)) to shorthand (e.g., H5N4A1).
    The order parameter allows specifying the order of the shorthand format."""
    # Remove any numbers after the % sign
    longhand = re.sub(r'%\s*\d+(\.\d+)?', '', longhand)
    matches = re.findall(r'(HexNAc|Hex|Fuc|NeuAc|NeuGc|Pent|Sulpho|Phospho|dHex)\((\d+)\)', longhand)
    default_order = ['HexNAc', 'Hex', 'Fuc', 'NeuAc', 'NeuGc', 'Pent', 'Phospho', 'Sulpho']
    # Replace 'dHex' with 'Fuc' in matches
    matches = [(name if name != 'dHex' else 'Fuc', count) for name, count in matches]
    order = order or default_order
    ordered_matches = sorted(matches, key=lambda x: order.index(x[0]) if x[0] in order else len(order))
    shorthand = ''.join(f"{LONGHAND_TO_SHORTHAND[name]}{count}" for name, count in ordered_matches)
    return shorthand
    
def detect_format(glycan):
    """Detects the format of a glycan string (shorthand or longhand)."""
    if re.match(r"^H\d+N\d+(A\d+|F\d+|G\d+|X\d+|A\d+)*$", glycan):  # Adjust regex for NeuGc, Xyl
        return "shorthand"
    elif re.match(r'^(HexNAc|Hex|Fuc|NeuAc|NeuGc|Pent|Sulpho|Phospho|dHex)\(\d+\)', glycan):
        return 'longhand'
    else:
        raise ValueError("Unknown glycan format.")
    raise ValueError("Unknown glycan format.")

def process_csv(input_file, output_file, glycan_column):
    """Reads a CSV file, detects the glycan format, and adds a new column for the converted glycans."""
    try:
        print(f"Reading file: {input_file}")
        df = pd.read_csv(input_file, dtype=str)  # Force all columns to be read as strings
        print(df.head())

        # Detect format of the first glycan
        first_glycan = df[glycan_column].iloc[0]
        print(f"First glycan detected: {first_glycan}")

        format_type = detect_format(first_glycan)
        print(f"Detected format: {format_type}")

        # Add a new column for the converted glycans
        if format_type == "shorthand":
            df['converted_glycan'] = df[glycan_column].apply(convert_shorthand_to_longhand)
        else:
            df['converted_glycan'] = df[glycan_column].apply(convert_longhand_to_shorthand)

        # Write the updated dataframe to the output file
        df.to_csv(output_file, index=False)

    except Exception as e:
        print("Error:", e)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert glycan formats between shorthand and longhand notation.")
    parser.add_argument("-i", "--input_file", required=True, help="Path to input CSV file.")
    parser.add_argument("-o", "--output_file", help="Path to output CSV file. Default is {input_filename}_converted.csv in the same directory as the input file.")
    parser.add_argument("-g", "--glycan_column", default="Glycans", help="Name of the glycan column in the CSV file (default: Glycans).")
    parser.add_argument("-d", "--direction", choices=["to_longhand", "to_shorthand"], help="Conversion direction: to_longhand (H5N4 -> HexNAc(4)Hex(5)) or to_shorthand (HexNAc(4)Hex(5) -> H5N4). If not provided, the format will be detected automatically.")
    args = parser.parse_args()

    # If no output file is provided, create a default output file name
    if not args.output_file:
        input_filename = os.path.splitext(os.path.basename(args.input_file))[0]
        input_dir = os.path.dirname(args.input_file)
        args.output_file = os.path.join(input_dir, f"{input_filename}_converted.csv")
    
    # Process the CSV file
    process_csv(args.input_file, args.output_file, args.glycan_column)
    print(f"File written to: {args.output_file}")
