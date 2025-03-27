# Glycopeptide Sequence Finder

A tool for identifying and analyzing glycopeptide sequences from protein FASTA files. 

Check out the [wiki](https://github.com/RichardDShipman/Glycopeptide_Sequence_Finder/wiki)!

# Quick Start Guide

## Project Structure

The project follows a standard Python package structure:

```
Glycopeptide_Sequence_Finder/
├── data/                           # Data directory containing all input/output files
│   ├── digested_peptide_library/  # Output peptide library
│   ├── digested_glycopeptide_library/  # Output glycopeptide library
│   ├── glycan_mass_library/       # Glycan mass library files
│   ├── glycosite_library/         # Output glycosite library
│   ├── logs/                      # Log files
│   ├── mock_mass_spectra/         # Mock mass spectrometry data
│   ├── glycan_library/            # Glycan library files
│   └── test_proteomes/            # Input FASTA files
├── src/                           # Source code directory
│   └── glycopeptide_sequence_finder_cmd.py  # Main script
├── scripts/                       # Utility scripts
│   └── batch_glycopeptide_sequence_finder.sh  # Batch processing script
├── Dockerfile                     # Container configuration
└── requirements.txt               # Python dependencies
```

## Installation

1. Clone the repository:
```bash
git clone https://github.com/RichardDShipman/Glycopeptide_Sequence_Finder.git
cd Glycopeptide_Sequence_Finder
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Command Line Interface

Run the script from the project root directory:

```bash
python src/glycopeptide_sequence_finder_cmd.py -i data/test_proteomes/your_input.fasta -o data/digested_glycopeptide_library/output.csv
```

### Batch Processing

Run the batch script to process multiple FASTA files:

```bash
./scripts/batch_glycopeptide_sequence_finder.sh
```

### Docker

Build the Docker image:
```bash
docker build -t glycopeptide-sequence-finder .
```

Run the container:
```bash
docker run -v "$(pwd)/data/test_proteomes:/app/data/test_proteomes" \
           -v "$(pwd)/data/digested_glycopeptide_library:/app/data/digested_glycopeptide_library" \
           glycopeptide-finder -i /app/data/test_proteomes/phosvitin_uniprotkb_2025_03_20.fasta \
           -o /app/data/digested_glycopeptide_library/phosvitin_output.csv -v
```

## Overview

Welcome to the Glycopeptide Sequence Finder! 

This tool aims to facilitate the theoretical study of glycopeptide sequences found in various species across the tree of life. It can be used with publicly available curated data from protein sequence databases or by uploading your own protein FASTA sequences files.

![MockMassSpectrumOfAN-GlycopeptideFromAPig!](/data/mock_mass_spectra/sp_P01042_KNG1_HUMAN_205_ITYSIVQTNCSK_G62765YT_mock_mass_spectrum.png)

*Example above*: Calculated mass spectrum for glycopeptide Kininogen-1 (P01042) KNG1_HUMAN - ITYSIVQTNCSK - 205 - G62765YT - HexNAc(2)Hex(8) created `plot_mock_mass_spectra.py`. Note: This is a basic N-glycan structure fragment under are minimum sequential fragmentation of Hex(8) down to HexHAc(2) core N-glycan.

**Glycopeptide Sequence Finder** is a Python script that processes protein sequences from a FASTA file to find amino acid sequences which may or may not contain the post-translational modification glycosylation, the attachment of glycans (polysaccharides) to protein sequences. It uses user-specified proteases to digest and cleave protein sequences into amino acid sequences. The script then identifies N-linked glycopeptides using glycosylation sequon (motifs) like the N-sequon "N[^P][STC]" (NX[STC], where X is not P), O-sequon "[S/T"], or C-sequon "W..[WCF]". It calculates the properties of these glycopeptides, including mass, hydrophobicity, and glycosylation sites. Additionally, the script gathers information from the inputted FASTA file to create a predicted digested glycopeptide (peptide sequence backbone) library. The output is written to a CSV file, making it easy to integrate into downstream analyses.

Additionally, a directory containing proteomics related info is stored in a directory titled `digested_peptide_library` for other use cases. This may be expanded in the future for other needs. Run the tool to generate, files too big for GitHub.

I am currently developing a basic calculation tool for fragment ions of high mannose N-glycans to simulate the fragmentation behavior of N-glycopeptides under HCD (High-Energy Collisional Dissociation) conditions. That is the image used above of a mock mass spectrum. This tool is designed to provide a fundamental outline of N-glycopeptide fragmentation, focusing on the generation of fragment ions typically observed in mass spectrometry experiments.

As part of the project, I am also integrating a plotting utility that visualizes the calculated mock fragment ion values, allowing for the comparison of the theoretical results against the mass spectrum data. This enables a clearer understanding of how glycopeptides and their glycans fragment during analysis.

For ease of experimentation and reproducibility, all data used in the calculations and plots, including example glycopeptides and their corresponding glycan compositions, are provided in organized folders. This setup allows for quick reference and testing of different theoretical models and fragmentation pathways.

*Note*: This project is for fun and more of an exploration of glycoproteomic space in silico. Who knows where this may lead.

# Table of Contents

[Overview](#overview)

[Requirements](#requirements)

[Installation](#installation)

[Usage](#usage)

### Reference Materials

[License](#license)

[Acknowledgments](#acknowledgments)

[Appendix](#appendix)

- [Batch run](#batch-processing-scripts)
- [Log File Details](#log-file)
- [Test Proteomes List](#test-proteomes)
- [Glycan Mass Library](#glycan-mass-library)
- [Test Proteomes](#test-proteomes)

## Features

1. **Protease-Specific Cleavage:**
    - Supports several commonly used proteases, including:
        - **Trypsin:** Cleaves after K or R, except if followed by P.
        - **Chymotrypsin:** Cleaves after F, W, or Y, except if followed by P.
        - **Glu-C:** Cleaves after E.
        - **Lys-C:** Cleaves after K.
        - **Arg-C:** Cleaves after R.
        - **Asp-N:** Cleaves before D.
        - **Pepsin:** Cleaves after F, L, W, or Y.
        - **Proteinase K:** Cleaves after A, F, I, L, V, W, or Y.
        - **All:** Runs all proteases above.
2. **Missed Cleavages:**
    - Allows specifying the number of missed cleavages to simulate incomplete digestion.
3. **Peptide Length Control:**
    - Minimum peptide length filter (default: 5 amino acids)
    - Maximum peptide length filter (default: 50 amino acids)
4. **Glycosylation Type:** 
    - Select from N-linked (N), O-linked (O), or C-linked (C) glycopeptides. (Adjust or add sequon)
        - **N-linked:** N-sequon "N[^P][STC]"
        - **O-linked:** O-sequon "[ST]" 
        - **C-linked:** C-sequon "W..[WCF]"
        - Others in the works.
5. **Peptide Property Calculation:**
    - Calculates peptide mass, hydrophobicity, isoelectric point (pI), charge states m/z values, and N-glycan ion fragmentation series (experimental).

## Requirements

- Python 3.7 or later
- Libraries:
    - argparse
    - csv
    - re
    - biopython

## Installation

Install the required Python libraries using pip:

```sh
pip install argparse biopython pandas
```

## Usage

Run the script from the command line with the following arguments:

```sh
python glycopeptide_finder_cmd.py -i <input_fasta> [-o <output_csv>] [-p <protease>] [-g <glycosylation>] [-c <missed_cleavages>] [-l <log.txt>] [-v]
```

### Arguments

- `-i`, `--input` (required): Path to the input FASTA file.
- `-o`, `--output` (optional): Path to the output CSV file. If omitted, a default name is generated.
- `-p`, `--protease` (optional): Protease to use for cleavage. Default is trypsin.
- `-g`, `--glycosylation` (optional): Glycosylation sequon to find in peptides. Default is N-linked. (N, O, C) Warning when using O or C, experimental.
- `-c`, `--missed_cleavages` (optional): Number of missed cleavages allowed. Default is 0.
- `-l log.txt`, `--log log.txt` (optional): Path to the log file. If omitted, logging is disabled.
- `-v`, `--verbose` (optional): Enable verbose output. Default is False.
- `-y`, `--glycan`: Path to the glycan file (CSV format) (Default, 4 glycans stored in file). 
- `-z`, `--charge`: (Optional) Maximum charge state to compute (default: 5).
- `-m`, `--max_peptide_length`: (Optional) Max peptide length after digestion (default: 50).
- `--min-length`: (Optional) Minimum peptide length after digestion (default: 5).

### Example

```sh
python glycopeptide_finder_cmd.py -i test_proteomes/human_uniprotkb_proteome_UP000005640_AND_revi_2025_01_17.fasta -p trypsin -g N -c 0 -z 3 --min-length 6
```

The output file will be dynamically named:

`example_predicted_trypsin_glycopeptides.csv`

### Example CSV Content

Ion series data can be generated with `--ion-series` arguement flag, example below. Prototype at the moment. Work in progress. Moving to separate function.

```CSV
ProteinID,Site,GlyToucan_AC,Composition,ShorthandGlycan,Peptide,Start,End,Length,Sequon,GlycopeptideMass,PeptideMass,GlycanMass,Hydrophobicity,pI,z2,Charge,IonSeries
sp|O95445|APOM_HUMAN,135.0,G22768VO,HexNAc(2)Hex(3),N2H3,TELFSSSCPGGIMLNETGQGYQR,121.0,143.0,23.0,NET,3690.543457999999,2474.1205949999994,1216.422863,-0.47826,4.26,1846.2790049999996,2,"{'b': [102.055, 231.0975, 344.1816, 491.25, 578.282, 665.3141, 752.3461, 855.3553, 952.4081, 1009.4295, 1066.451, 1179.535, 1310.5755, 1423.6596, 1537.7025, 1666.7451, 1767.7928, 1824.8142, 1952.8728, 2009.8943, 2172.9576, 2301.0162], 'y': [175.119, 303.1775, 466.2409, 523.2623, 651.3209, 708.3424, 809.39, 938.4326, 1052.4756, 1165.5596, 1296.6001, 1409.6842, 1466.7056, 1523.7271, 1620.7799, 1723.789, 1810.8211, 1897.8531, 1984.8851, 2131.9535, 2245.0376, 2374.0802], 'c': [119.0815, 248.124, 361.2081, 508.2765, 595.3085, 682.3406, 769.3726, 872.3818, 969.4346, 1026.456, 1083.4775, 1196.5615, 1327.602, 1440.6861, 1554.729, 1683.7716, 1784.8193, 1841.8407, 1969.8993, 2026.9208, 2189.9841, 2318.0427], 'z': [140.0819, 268.1405, 431.2038, 488.2253, 616.2838, 673.3053, 774.353, 903.3956, 1017.4385, 1130.5226, 1261.563, 1374.6471, 1431.6686, 1488.69, 1585.7428, 1688.752, 1775.784, 1862.816, 1949.8481, 2096.9165, 2210.0005, 2339.0431], 'Y': {'Y0': 2475.1279, 'Y1': 2678.2073, 'Y2': 2881.2867, 'Y3': 3043.3395, 'Y4': 3205.3923, 'Y5': 3367.4451}, '2Y': {'2Y0': 1237.5639, '2Y1': 1339.1036, '2Y2': 1440.6433, '2Y3': 1521.6697, '2Y4': 1602.6961, '2Y5': 1683.7225}, 'B': {'B_HexNAc_1': 204.0867, 'B_HexNAc_2': 407.1661, 'B_Hex_1': 569.2189, 'B_Hex_2': 731.2717, 'B_Hex_3': 893.3245}, 'oxonium': {'ox_HexNAc': 204.0867, 'ox_Hex': 163.0601}}"
sp|P00450|CERU_HUMAN,138.0,G22768VO,HexNAc(2)Hex(3),N2H3,EHEGAIYPDNTTDFQR,129.0,144.0,16.0,NTT,3108.2565080000004,1891.8336450000002,1216.422863,-1.51875,3.95,1555.1355300000002,2,"{'b': [130.0499, 267.1088, 396.1514, 453.1728, 524.2099, 637.294, 800.3573, 897.4101, 1012.437, 1126.48, 1227.5276, 1328.5753, 1443.6023, 1590.6707, 1718.7292], 'y': [175.119, 303.1775, 450.2459, 565.2729, 666.3206, 767.3682, 881.4112, 996.4381, 1093.4909, 1256.5542, 1369.6383, 1440.6754, 1497.6968, 1626.7394, 1763.7983], 'c': [147.0764, 284.1353, 413.1779, 470.1993, 541.2364, 654.3205, 817.3838, 914.4366, 1029.4635, 1143.5065, 1244.5541, 1345.6018, 1460.6288, 1607.6972, 1735.7557], 'z': [140.0819, 268.1405, 415.2089, 530.2358, 631.2835, 732.3312, 846.3741, 961.401, 1058.4538, 1221.5171, 1334.6012, 1405.6383, 1462.6598, 1591.7024, 1728.7613], 'Y': {'Y0': 1892.8409, 'Y1': 2095.9203, 'Y2': 2298.9997, 'Y3': 2461.0525, 'Y4': 2623.1053, 'Y5': 2785.1581}, '2Y': {'2Y0': 946.4205, '2Y1': 1047.9602, '2Y2': 1149.4999, '2Y3': 1230.5263, '2Y4': 1311.5527, '2Y5': 1392.5791}, 'B': {'B_HexNAc_1': 204.0867, 'B_HexNAc_2': 407.1661, 'B_Hex_1': 569.2189, 'B_Hex_2': 731.2717, 'B_Hex_3': 893.3245}, 'oxonium': {'ox_HexNAc': 204.0867, 'ox_Hex': 163.0601}}"
sp|P00450|CERU_HUMAN,227.0,G22768VO,HexNAc(2)Hex(3),N2H3,EFVVMFSVVDENFSWYLEDNIK,216.0,237.0,22.0,NFS,3925.6900780000005,2709.2672150000003,1216.422863,0.14545,3.39,1963.8523150000003,2,"{'b': [130.0499, 277.1183, 376.1867, 475.2551, 606.2956, 753.364, 840.396, 939.4644, 1038.5328, 1153.5598, 1282.6024, 1396.6453, 1543.7137, 1630.7457, 1816.8251, 1979.8884, 2092.9724, 2222.015, 2337.042, 2451.0849, 2564.169], 'y': [147.1128, 260.1969, 374.2398, 489.2667, 618.3093, 731.3934, 894.4567, 1080.536, 1167.5681, 1314.6365, 1428.6794, 1557.722, 1672.7489, 1771.8173, 1870.8857, 1957.9178, 2104.9862, 2236.0267, 2335.0951, 2434.1635, 2581.2319], 'c': [147.0764, 294.1448, 393.2132, 492.2816, 623.3221, 770.3905, 857.4225, 956.4909, 1055.5593, 1170.5863, 1299.6289, 1413.6718, 1560.7402, 1647.7722, 1833.8516, 1996.9149, 2109.9989, 2239.0415, 2354.0685, 2468.1114, 2581.1955], 'z': [112.0757, 225.1598, 339.2027, 454.2297, 583.2723, 696.3563, 859.4196, 1045.499, 1132.531, 1279.5994, 1393.6423, 1522.6849, 1637.7119, 1736.7803, 1835.8487, 1922.8807, 2069.9491, 2200.9896, 2300.058, 2399.1264, 2546.1948], 'Y': {'Y0': 2710.2745, 'Y1': 2913.3539, 'Y2': 3116.4333, 'Y3': 3278.4861, 'Y4': 3440.5389, 'Y5': 3602.5917}, '2Y': {'2Y0': 1355.1372, '2Y1': 1456.6769, '2Y2': 1558.2166, '2Y3': 1639.243, '2Y4': 1720.2694, '2Y5': 1801.2958}, 'B': {'B_HexNAc_1': 204.0867, 'B_HexNAc_2': 407.1661, 'B_Hex_1': 569.2189, 'B_Hex_2': 731.2717, 'B_Hex_3': 893.3245}, 'oxonium': {'ox_HexNAc': 204.0867, 'ox_Hex': 163.0601}}"
```

## Protease Rules

The following proteases are supported:

| Protease      | Cleavage Rule                        |
|---------------|--------------------------------------|
| Trypsin       | After K or R, not P                  |
| Chymotrypsin  | After F, W, or Y, not P              |
| Glu-C         | After E                              |
| Lys-C         | After K                              |
| Arg-C         | After R                              |
| Asp-N         | Before D                             |
| Pepsin        | After F, L, W, or Y                  |
| Proteinase K  | After A, F, I, L, V, W, or Y         |
| All           | Runs all proteases above             |

## Glycosylation Type Rules

The following glycosylation types sequons (motifs) are supported:

| Glycosylation Type | Sequon Pattern |
|--------------------|----------------|
| N-linked           | N[^P][STC]     |
| O-linked           | [ST]           |
| C-linked           | W..[WCF]       |

## Glycan Library

The default glycan mass library is defined as a DataFrame containing a set of glycans with their respective compositions and masses. This library is used to calculate the properties of glycopeptides. Alter if you wish to change the glycan mass library in the script

```python
default_glycan_library = pd.DataFrame([
    #{"glytoucan_ac": "G80920RR", "byonic": "HexNAc(2)Hex(9) % 1864.634157", "composition": "HexNAc(2)Hex(9)", "mass": 1864.634157}, # N2H9
    {"glytoucan_ac": "G62765YT", "byonic": "HexNAc(2)Hex(8) % 1702.581333", "composition": "HexNAc(2)Hex(8)", "mass": 1702.581333}, # N2H8
    #{"glytoucan_ac": "G31852PQ", "byonic": "HexNAc(2)Hex(7) % 1540.528510", "composition": "HexNAc(2)Hex(7)", "mass": 1540.528510}, # N2H7
    #{"glytoucan_ac": "G41247ZX", "byonic": "HexNAc(2)Hex(6) % 1378.475686", "composition": "HexNAc(2)Hex(6)", "mass": 1378.475686}, # N2H6
])
```

Other libraries in file. Feel free to expand to meet needs of user.

This DataFrame includes the following columns:
- `glytoucan_ac`: The glycosylation context identifier.
- `byonic`: The peptide sequence and mass in Byonic format.
- `composition`: The glycan composition.
- `mass`: The mass of the glycan.

The default glycan mass library can be expanded or customized as needed for specific analyses.

### Example glycan mass library data

This data is stored in the `glycan_mass_library` directory.

```csv
glytoucan_ac,byonic,composition,converted_glycan,glycan_composition_sequence
G62765YT,HexNAc(2)Hex(8) % 1702.581333,HexNAc(2)Hex(8),N2H8,NNHHHHHHHH
G31852PQ,HexNAc(2)Hex(7) % 1540.528510,HexNAc(2)Hex(7),N2H7,NNHHHHHHH
G41247ZX,HexNAc(2)Hex(6) % 1378.475686,HexNAc(2)Hex(6),N2H6,NNHHHHHH
```

# Plot Mock Mass Spectrum

![PiggieGlycopeptides](/data/mock_mass_spectra/sp_P02763_A1AG1_HUMAN_103_ENGTISR_G62765YT_mock_mass_spectrum.png)

plot_mock_mass_spectra.py

### *Usage Guide*

*This script generates mock mass spectrum plots from glycopeptide ion series stored in a CSV file. It processes ion series data, assigns colors to different ion types, computes ion numbers, and generates labeled mass spectrum plots.*

- Ion labels are automatically assigned based on their types.
- If a peptide sequence exceeds 50 characters, it will be skipped.
- The script ensures unique filenames for each output image.

#### Arguments:

- `-i, --input` (required): Path to the input CSV file containing glycopeptide ion series data.
- `-o, --output` (optional): Directory to save the generated plots (default: `mock_mass_spectra` directory).

### Input CSV Format

The CSV file must contain the following required columns:

- `IonSeries`: Dictionary-like string containing ion data (b, y, Y, B, oxonium).
- `ProteinID`: Identifier for the protein.
- `Peptide`: Peptide sequence.
- `Composition`: Glycan composition.
- `GlyToucan_AC`: GlyToucan accession number.

```bash
python script.py -i digested_glycopeptide_library/pig_uniprotkb_proteome_UP000008227_AND_revi_2025_02_01_trypsin_digested_mc0_z2_N-glycopeptides.csv
```

This will process `example_data.csv` and save the mass spectrum plots in the `results/` directory.

## Glycan Composition to Sequence Converter

This Python script converts glycan composition data into a glycan sequence. It takes an input CSV file containing a column with glycan composition (e.g., N2H3F1) and generates a glycan sequence, expanding the monosaccharides according to their counts. The results are saved to a new CSV file.

Usage

```bash
python glycan_converter.py -i <input_file> -o <output_file> -g <glycan_column>
```
- -i, --input: Path to the input CSV file.
- -o, --output: Path to save the output CSV file (default: input filename with "_composition_sequences.csv" suffix).
- -g, --glycan: The name of the column containing glycan composition data.

### Requirements
	•	pandas
	•	re (built-in Python module)

### Example

```bash
python glycan_converter.py -i glycans.csv -g glycan_composition
```

## Create Glycan Mass Library 

create_glycan_library.py

This script processes glycan data from a CSV file and splits it into columns based on specific formatting rules.

## Overview
The script takes an input CSV file with glyc肽 data formatted in two columns:
1. **glytoucan_ac**: A string representing the glycosylation context.
2. **byonic** / **byonic_sequence**: A composite column containing two pieces of information separated by a '%' character:
    - The peptide sequence (composition).
    - A numerical value representing mass.

The script splits the `byonic` column into its constituent parts, creating three new columns in the output file:
1. **glytoucan_ac**: The glycosylation context.
2. **composition**: The peptide sequence.
3. **mass**: The numerical mass value.

## Key Features

### Input File Format
The input CSV file must have exactly two columns per row, with the second column formatted as `<peptide> % <mass>`.

Example:
```
"glytoucan_ac","sequence_byonic","name_source"
"G00002CF","Hex(2)NeuGc(2) % 956.29687423","G93218EI"
"G00009BX","HexNAc(2)Hex(2)dHex(1) % 894.33286578","G93579XB"
"G00012MO","HexNAc(1)Hex(3) % 707.248407805","G08590QR"
```

### Output File Format
The output file will have three columns:
- The first column (`glytoucan_ac`) remains unchanged.
- The second column contains the peptide sequence from `byonic`.
- The third column contains the numerical mass value.

Example:
```
glytoucan_ac,byonic,composition,converted_glycan,glycan_composition_sequence,composition,mass
G62765YT,HexNAc(2)Hex(8) % 1702.581333,HexNAc(2)Hex(8),N2H8,NNHHHHHHHH,HexNAc(2)Hex(8),1702.581333
G31852PQ,HexNAc(2)Hex(7) % 1540.528510,HexNAc(2)Hex(7),N2H7,NNHHHHHHH,HexNAc(2)Hex(7),1540.528510
G41247ZX,HexNAc(2)Hex(6) % 1378.475686,HexNAc(2)Hex(6),N2H6,NNHHHHHH,HexNAc(2)Hex(6),1378.475686
```

The output CSV will contain glycans ranked by their weighted adjusted HF score.

# Batch Processing Scripts

Shell scripts for batch processing.

## Batch Run for FASTA Processing

batch_glycopeptide_sequence_finder.sh

To process multiple FASTA files in parallel using all proteases, run the following command:

```sh
./batch_glycopeptide_sequence_finder.sh
```

Parameters can be adjusted in the shell script.

### Parameters

- `ls test_proteomes/*.fasta`: Lists all FASTA files in the `test_proteomes` directory.
- `xargs -I {} -P 4`: Executes the command in parallel with up to 4 processes. The `{}` is a placeholder for each file name.
- `python glycopeptide_finder_cmd.py`: The script to run for each FASTA file.
- `-i "{}"`: Specifies the input FASTA file, where `{}` is replaced by each file name.
- `-p all`: Uses all proteases for cleavage.
- `-g N`: Searches for N-linked glycosylation sequons.
- `-c 0`: Allows 0 missed cleavages.
- `-v`: Enables verbose output.

This command allows you to efficiently process multiple FASTA files in parallel, reducing the overall processing time.

## Merging CSV Files

The script includes a function to merge all CSV files from a specified directory into a single CSV file. This can be useful for consolidating the results of multiple digestions into one file for easier analysis.

```sh
python merge_digested_glycopeptide_library.py
```

# Machine Learning (experimental)

Explore machine learning space with the glycopeptide data above. Experimental, in development.

## Glycopeptide One-Hot Encoding Script

This script encodes glycopeptide data into one-hot encoded feature vectors. It handles peptide sequences, glycan compositions, and charge states, generating a feature vector for each glycopeptide in the input data. This will be utilized for future Y class labels in glycoproteomics machine learning.

### Features
- Peptide Encoding: Encodes the 20 standard amino acids into a 50x20 matrix.
- Glycan Encoding: Encodes a simplified set of monosaccharides (N, H, F, A) into a 30x4 matrix.
- Charge State Encoding: Encodes charge states (1-10) into a one-hot vector.
- CSV Output: Generates an encoded dataset and an encoding definition CSV file.

## Usage

```bash
python encode_glycopeptides.py -i input.csv -o output.csv -d encoding_definition.csv
```

- -i : Input CSV file with peptide, glycan, and charge data.
- -o : Output CSV file for encoded data.
- -d : Output CSV for encoding definitions.

### Dependencies
- numpy
- pandas
- argparse
- re

### Example

```sh
python glycopeptide_finder_cmd.py -i test_proteomes/human_uniprotkb_proteome_UP000005640_AND_revi_2025_01_17.fasta -p trypsin -g N -c 0 -z 3
```

The output file will be dynamically named:

`example_predicted_trypsin_glycopeptides.csv`

### Example CSV Content

```CSV
ProteinID,Site,GlyToucan_AC,Composition,ShorthandGlycan,Peptide,Start,End,Length,Sequon,GlycopeptideMass,PeptideMass,GlycanMass,Hydrophobicity,pI,z2,Charge,IonSeries,Glycan_Composition_Sequence,One_Hot_Encoding
ProteinID,Site,GlyToucan_AC,Composition,ShorthandGlycan,Peptide,Start,End,Length,Sequon,GlycopeptideMass,PeptideMass,GlycanMass,Hydrophobicity,pI,z2,Charge,IonSeries
sp|O95445|APOM_HUMAN,135.0,G22768VO,HexNAc(2)Hex(3),N2H3,TELFSSSCPGGIMLNETGQGYQR,121.0,143.0,23.0,NET,3690.543457999999,2474.1205949999994,1216.422863,-0.47826,4.26,1846.2790049999996,2,"{'b': [102.055, 231.0975, 344.1816, 491.25, 578.282, 665.3141, 752.3461, 855.3553, 952.4081, 1009.4295, 1066.451, 1179.535, 1310.5755, 1423.6596, 1537.7025, 1666.7451, 1767.7928, 1824.8142, 1952.8728, 2009.8943, 2172.9576, 2301.0162], 'y': [175.119, 303.1775, 466.2409, 523.2623, 651.3209, 708.3424, 809.39, 938.4326, 1052.4756, 1165.5596, 1296.6001, 1409.6842, 1466.7056, 1523.7271, 1620.7799, 1723.789, 1810.8211, 1897.8531, 1984.8851, 2131.9535, 2245.0376, 2374.0802], 'c': [119.0815, 248.124, 361.2081, 508.2765, 595.3085, 682.3406, 769.3726, 872.3818, 969.4346, 1026.456, 1083.4775, 1196.5615, 1327.602, 1440.6861, 1554.729, 1683.7716, 1784.8193, 1841.8407, 1969.8993, 2026.9208, 2189.9841, 2318.0427], 'z': [140.0819, 268.1405, 431.2038, 488.2253, 616.2838, 673.3053, 774.353, 903.3956, 1017.4385, 1130.5226, 1261.563, 1374.6471, 1431.6686, 1488.69, 1585.7428, 1688.752, 1775.784, 1862.816, 1949.8481, 2096.9165, 2210.0005, 2339.0431], 'Y': {'Y0': 2475.1279, 'Y1': 2678.2073, 'Y2': 2881.2867, 'Y3': 3043.3395, 'Y4': 3205.3923, 'Y5': 3367.4451}, '2Y': {'2Y0': 1237.5639, '2Y1': 1339.1036, '2Y2': 1440.6433, '2Y3': 1521.6697, '2Y4': 1602.6961, '2Y5': 1683.7225}, 'B': {'B_HexNAc_1': 204.0867, 'B_HexNAc_2': 407.1661, 'B_Hex_1': 569.2189, 'B_Hex_2': 731.2717, 'B_Hex_3': 893.3245}, 'oxonium': {'ox_HexNAc': 204.0867, 'ox_Hex': 163.0601}}"
sp|P00450|CERU_HUMAN,138.0,G22768VO,HexNAc(2)Hex(3),N2H3,EHEGAIYPDNTTDFQR,129.0,144.0,16.0,NTT,3108.2565080000004,1891.8336450000002,1216.422863,-1.51875,3.95,1555.1355300000002,2,"{'b': [130.0499, 267.1088, 396.1514, 453.1728, 524.2099, 637.294, 800.3573, 897.4101, 1012.437, 1126.48, 1227.5276, 1328.5753, 1443.6023, 1590.6707, 1718.7292], 'y': [175.119, 303.1775, 450.2459, 565.2729, 666.3206, 767.3682, 881.4112, 996.4381, 1093.4909, 1256.5542, 1369.6383, 1440.6754, 1497.6968, 1626.7394, 1763.7983], 'c': [147.0764, 284.1353, 413.1779, 470.1993, 541.2364, 654.3205, 817.3838, 914.4366, 1029.4635, 1143.5065, 1244.5541, 1345.6018, 1460.6288, 1607.6972, 1735.7557], 'z': [140.0819, 268.1405, 415.2089, 530.2358, 631.2835, 732.3312, 846.3741, 961.401, 1058.4538, 1221.5171, 1334.6012, 1405.6383, 1462.6598, 1591.7024, 1728.7613], 'Y': {'Y0': 1892.8409, 'Y1': 2095.9203, 'Y2': 2298.9997, 'Y3': 2461.0525, 'Y4': 2623.1053, 'Y5': 2785.1581}, '2Y': {'2Y0': 946.4205, '2Y1': 1047.9602, '2Y2': 1149.4999, '2Y3': 1230.5263, '2Y4': 1311.5527, '2Y5': 1392.5791}, 'B': {'B_HexNAc_1': 204.0867, 'B_HexNAc_2': 407.1661, 'B_Hex_1': 569.2189, 'B_Hex_2': 731.2717, 'B_Hex_3': 893.3245}, 'oxonium': {'ox_HexNAc': 204.0867, 'ox_Hex': 163.0601}}"
sp|P00450|CERU_HUMAN,227.0,G22768VO,HexNAc(2)Hex(3),N2H3,EFVVMFSVVDENFSWYLEDNIK,216.0,237.0,22.0,NFS,3925.6900780000005,2709.2672150000003,1216.422863,0.14545,3.39,1963.8523150000003,2,"{'b': [130.0499, 277.1183, 376.1867, 475.2551, 606.2956, 753.364, 840.396, 939.4644, 1038.5328, 1153.5598, 1282.6024, 1396.6453, 1543.7137, 1630.7457, 1816.8251, 1979.8884, 2092.9724, 2222.015, 2337.042, 2451.0849, 2564.169], 'y': [147.1128, 260.1969, 374.2398, 489.2667, 618.3093, 731.3934, 894.4567, 1080.536, 1167.5681, 1314.6365, 1428.6794, 1557.722, 1672.7489, 1771.8173, 1870.8857, 1957.9178, 2104.9862, 2236.0267, 2335.0951, 2434.1635, 2581.2319], 'c': [147.0764, 294.1448, 393.2132, 492.2816, 623.3221, 770.3905, 857.4225, 956.4909, 1055.5593, 1170.5863, 1299.6289, 1413.6718, 1560.7402, 1647.7722, 1833.8516, 1996.9149, 2109.9989, 2239.0415, 2354.0685, 2468.1114, 2581.1955], 'z': [112.0757, 225.1598, 339.2027, 454.2297, 583.2723, 696.3563, 859.4196, 1045.499, 1132.531, 1279.5994, 1393.6423, 1522.6849, 1637.7119, 1736.7803, 1835.8487, 1922.8807, 2069.9491, 2200.9896, 2300.058, 2399.1264, 2546.1948], 'Y': {'Y0': 2710.2745, 'Y1': 2913.3539, 'Y2': 3116.4333, 'Y3': 3278.4861, 'Y4': 3440.5389, 'Y5': 3602.5917}, '2Y': {'2Y0': 1355.1372, '2Y1': 1456.6769, '2Y2': 1558.2166, '2Y3': 1639.243, '2Y4': 1720.2694, '2Y5': 1801.2958}, 'B': {'B_HexNAc_1': 204.0867, 'B_HexNAc_2': 407.1661, 'B_Hex_1': 569.2189, 'B_Hex_2': 731.2717, 'B_Hex_3': 893.3245}, 'oxonium': {'ox_HexNAc': 204.0867, 'ox_Hex': 163.0601}}"
```

## Glycan Library

The default glycan mass library is defined as a DataFrame containing a set of glycans with their respective compositions and masses. This library is used to calculate the properties of glycopeptides. Alter if you wish to change the glycan mass library in the script

```python
default_glycan_library = pd.DataFrame([
    #{"glytoucan_ac": "G80920RR", "byonic": "HexNAc(2)Hex(9) % 1864.634157", "composition": "HexNAc(2)Hex(9)", "mass": 1864.634157}, # N2H9
    {"glytoucan_ac": "G62765YT", "byonic": "HexNAc(2)Hex(8) % 1702.581333", "composition": "HexNAc(2)Hex(8)", "mass": 1702.581333}, # N2H8
    #{"glytoucan_ac": "G31852PQ", "byonic": "HexNAc(2)Hex(7) % 1540.528510", "composition": "HexNAc(2)Hex(7)", "mass": 1540.528510}, # N2H7
    #{"glytoucan_ac": "G41247ZX", "byonic": "HexNAc(2)Hex(6) % 1378.475686", "composition": "HexNAc(2)Hex(6)", "mass": 1378.475686}, # N2H6
])
```

Other libraries in file. Feel free to expand to meet needs of user.

This DataFrame includes the following columns:
- `glytoucan_ac`: The glycosylation context identifier.
- `byonic`: The peptide sequence and mass in Byonic format.
- `composition`: The glycan composition.
- `mass`: The mass of the glycan.

The default glycan mass library can be expanded or customized as needed for specific analyses.

### Example glycan mass library data

This data is stored in the `glycan_mass_library` directory.

```csv
glytoucan_ac,byonic,composition,converted_glycan,glycan_composition_sequence
G62765YT,HexNAc(2)Hex(8) % 1702.581333,HexNAc(2)Hex(8),N2H8,NNHHHHHHHH
G31852PQ,HexNAc(2)Hex(7) % 1540.528510,HexNAc(2)Hex(7),N2H7,NNHHHHHHH
G41247ZX,HexNAc(2)Hex(6) % 1378.475686,HexNAc(2)Hex(6),N2H6,NNHHHHHH
```

# Plot Mock Mass Spectrum

![PiggieGlycopeptides](/data/mock_mass_spectra/sp_P02763_A1AG1_HUMAN_103_ENGTISR_G62765YT_mock_mass_spectrum.png)

plot_mock_mass_spectra.py

### *Usage Guide*

*This script generates mock mass spectrum plots from glycopeptide ion series stored in a CSV file. It processes ion series data, assigns colors to different ion types, computes ion numbers, and generates labeled mass spectrum plots.*

- Ion labels are automatically assigned based on their types.
- If a peptide sequence exceeds 50 characters, it will be skipped.
- The script ensures unique filenames for each output image.

#### Arguments:

- `-i, --input` (required): Path to the input CSV file containing glycopeptide ion series data.
- `-o, --output` (optional): Directory to save the generated plots (default: `mock_mass_spectra` directory).

### Input CSV Format

The CSV file must contain the following required columns:

- `IonSeries`: Dictionary-like string containing ion data (b, y, Y, B, oxonium).
- `ProteinID`: Identifier for the protein.
- `Peptide`: Peptide sequence.
- `Composition`: Glycan composition.
- `GlyToucan_AC`: GlyToucan accession number.

```bash
python plot_mock_mass_spectra.py -i data/digested_glycopeptide_library/pig_uniprotkb_proteome_UP000008227_AND_revi_2025_02_01_trypsin_digested_mc0_z2_N-glycopeptides.csv
```

This will process `example_data.csv` and save the mass spectrum plots in the `results/` directory.

## Glycan Composition to Sequence Converter

This Python script converts glycan composition data into a glycan sequence. It takes an input CSV file containing a column with glycan composition (e.g., N2H3F1) and generates a glycan sequence, expanding the monosaccharides according to their counts. The results are saved to a new CSV file.

Usage

```bash
python glycan_converter.py -i <input_file> -o <output_file> -g <glycan_column>
```
- -i, --input: Path to the input CSV file.
- -o, --output: Path to save the output CSV file (default: input filename with “_composition_sequences.csv” suffix).
- -g, --glycan: The name of the column containing glycan composition data.

### Requirements
	•	pandas
	•	re (built-in Python module)

### Example

```bash
python glycan_converter.py -i glycans.csv -g glycan_composition
```

## Create Glycan Mass Library 

create_glycan_library.py

This script processes glycan data from a CSV file and splits it into columns based on specific formatting rules.

## Overview
The script takes an input CSV file with glyc肽 data formatted in two columns:
1. **glytoucan_ac**: A string representing the glycosylation context.
2. **byonic** / **byonic_sequence**: A composite column containing two pieces of information separated by a '%' character:
    - The peptide sequence (composition).
    - A numerical value representing mass.

The script splits the `byonic` column into its constituent parts, creating three new columns in the output file:
1. **glytoucan_ac**: The glycosylation context.
2. **composition**: The peptide sequence.
3. **mass**: The numerical mass value.

## Key Features

### Input File Format
The input CSV file must have exactly two columns per row, with the second column formatted as `<peptide> % <mass>`.

Example:
```
"glytoucan_ac","sequence_byonic","name_source"
"G00002CF","Hex(2)NeuGc(2) % 956.29687423","G93218EI"
"G00009BX","HexNAc(2)Hex(2)dHex(1) % 894.33286578","G93579XB"
"G00012MO","HexNAc(1)Hex(3) % 707.248407805","G08590QR"
```

### Output File Format
The output file will have three columns:
- The first column (`glytoucan_ac`) remains unchanged.
- The second column contains the peptide sequence from `byonic`.
- The third column contains the numerical mass value.

Example:
```
glytoucan_ac,byonic,composition,converted_glycan,glycan_composition_sequence,composition,mass
G62765YT,HexNAc(2)Hex(8) % 1702.581333,HexNAc(2)Hex(8),N2H8,NNHHHHHHHH,HexNAc(2)Hex(8),1702.581333
G31852PQ,HexNAc(2)Hex(7) % 1540.528510,HexNAc(2)Hex(7),N2H7,NNHHHHHHH,HexNAc(2)Hex(7),1540.528510
G41247ZX,HexNAc(2)Hex(6) % 1378.475686,HexNAc(2)Hex(6),N2H6,NNHHHHHH,HexNAc(2)Hex(6),1378.475686
```

The output CSV will contain glycans ranked by their weighted adjusted HF score.

# Batch Processing Scripts

Shell scripts for batch processing.

## Batch Run for FASTA Processing

batch_glycopeptide_sequence_finder.sh

To process multiple FASTA files in parallel using all proteases, run the following command:

```sh
./batch_glycopeptide_sequence_finder.sh
```

Parameters can be adjusted in the shell script.

### Parameters

- `ls test_proteomes/*.fasta`: Lists all FASTA files in the `test_proteomes` directory.
- `xargs -I {} -P 4`: Executes the command in parallel with up to 4 processes. The `{}` is a placeholder for each file name.
- `python glycopeptide_finder_cmd.py`: The script to run for each FASTA file.
- `-i "{}"`: Specifies the input FASTA file, where `{}` is replaced by each file name.
- `-p all`: Uses all proteases for cleavage.
- `-g N`: Searches for N-linked glycosylation sequons.
- `-c 0`: Allows 0 missed cleavages.
- `-v`: Enables verbose output.

This command allows you to efficiently process multiple FASTA files in parallel, reducing the overall processing time.

## Merging CSV Files

The script includes a function to merge all CSV files from a specified directory into a single CSV file. This can be useful for consolidating the results of multiple digestions into one file for easier analysis.

```sh
python merge_digested_glycopeptide_library.py
```

## Dockerfile

- Docker Setup for Glycopeptide Sequence Finder

This section explains how to build and run the Docker container for the Glycopeptide Sequence Finder.

1. Build the Docker Image

To create the Docker image, run the following command in the directory containing your Dockerfile and requirements.txt:

```sh
docker build -t gsf .
```

This will:
- Use the official Python 3.10-slim image.
- Set /app as the working directory.
- Install dependencies from requirements.txt.
- Copy the glycopeptide_sequence_finder_cmd.py script into the container.
- Set the entrypoint so that the script can be executed with arguments.

2. Run the Docker Container

To execute the script with test data, use:

```sh
docker run --rm \
    -v "$(pwd)/test_proteomes:/app/test_proteomes" \
    -v "$(pwd)/output:/app/digested_glycopeptide_library" \
    gsf \
    -i test_proteomes/apple_uniprotkb_proteome_UP000290289_AND_revi_2025_02_04.fasta \
    -g N \
    -o digested_glycopeptide_library/test.csv \
    -p chymotrypsin \
    -c 0 \
    -v
```

Explanation of Flags:
- --rm → Removes the container after execution.
- -v "$(pwd)/test_proteomes:/app/test_proteomes" → Mounts the input FASTA files.
- -v "$(pwd)/output:/app/digested_glycopeptide_library" → Mounts the output directory.
- gsf → Runs the built image.
- -i → Specifies the input FASTA file.
- -g → Sets the glycosylation type (default: N).
- -o → Defines the output file.
- -p → Specifies the protease (e.g., chymotrypsin).
- -c → Defines the missed cleavages.
- -v → Enables verbose mode.

3. Access the Output

The output files will be saved in the mounted directory on your local machine:

```sh
ls output/digested_glycopeptide_library/
```

Your results should be inside output/digested_glycopeptide_library/test.csv.

# Machine Learning (experimental)

Explore machine learning space with the glycopeptide data above. Experimental, in development.

## Glycopeptide One-Hot Encoding Script

This script encodes glycopeptide data into one-hot encoded feature vectors. It handles peptide sequences, glycan compositions, and charge states, generating a feature vector for each glycopeptide in the input data. This will be utilized for future Y class labels in glycoproteomics machine learning.

### Features
- Peptide Encoding: Encodes the 20 standard amino acids into a 50x20 matrix.
- Glycan Encoding: Encodes a simplified set of monosaccharides (N, H, F, A) into a 30x4 matrix.
- Charge State Encoding: Encodes charge states (1-10) into a one-hot vector.
- CSV Output: Generates an encoded dataset and an encoding definition CSV file.

## Usage

```bash
python encode_glycopeptides.py -i input.csv -o output.csv -d encoding_definition.csv
```

- -i : Input CSV file with peptide, glycan, and charge data.
- -o : Output CSV file for encoded data.
- -d : Output CSV for encoding definitions.

### Dependencies
- numpy
- pandas
- argparse
- re

## Output

```CSV
ProteinID,Site,GlyToucan_AC,Composition,ShorthandGlycan,Peptide,Start,End,Length,Sequon,GlycopeptideMass,PeptideMass,GlycanMass,Hydrophobicity,pI,z2,Charge,IonSeries,Glycan_Composition_Sequence,One_Hot_Encoding
sp|O95445|APOM_HUMAN,135.0,G22768VO,HexNAc(2)Hex(3),N2H3,TELFSSSCPGGIMLNETGQGYQR,121.0,143.0,23.0,NET,3690.543457999999,2474.1205949999994,1216.422863,-0.47826,4.26,1846.2790049999996,2,"{'b': [102.055, 231.0975, 344.1816, 491.25, 578.282, 665.3141, 752.3461, 855.3553, 952.4081, 1009.4295, 1066.451, 1179.535, 1310.5755, 1423.6596, 1537.7025, 1666.7451, 1767.7928, 1824.8142, 1952.8728, 2009.8943, 2172.9576, 2301.0162], 'y': [175.119, 303.1775, 466.2409, 523.2623, 651.3209, 708.3424, 809.39, 938.4326, 1052.4756, 1165.5596, 1296.6001, 1409.6842, 1466.7056, 1523.7271, 1620.7799, 1723.789, 1810.8211, 1897.8531, 1984.8851, 2131.9535, 2245.0376, 2374.0802], 'c': [119.0815, 248.124, 361.2081, 508.2765, 595.3085, 682.3406, 769.3726, 872.3818, 969.4346, 1026.456, 1083.4775, 1196.5615, 1327.602, 1440.6861, 1554.729, 1683.7716, 1784.8193, 1841.8407, 1969.8993, 2026.9208, 2189.9841, 2318.0427], 'z': [140.0819, 268.1405, 431.2038, 488.2253, 616.2838, 673.3053, 774.353, 903.3956, 1017.4385, 1130.5226, 1261.563, 1374.6471, 1431.6686, 1488.69, 1585.7428, 1688.752, 1775.784, 1862.816, 1949.8481, 2096.9165, 2210.0005, 2339.0431], 'Y': {'Y0': 2475.1279, 'Y1': 2678.2073, 'Y2': 2881.2867, 'Y3': 3043.3395, 'Y4': 3205.3923, 'Y5': 3367.4451}, '2Y': {'2Y0': 1237.5639, '2Y1': 1339.1036, '2Y2': 1440.6433, '2Y3': 1521.6697, '2Y4': 1602.6961, '2Y5': 1683.7225}, 'B': {'B_HexNAc_1': 204.0867, 'B_HexNAc_2': 407.1661, 'B_Hex_1': 569.2189, 'B_Hex_2': 731.2717, 'B_Hex_3': 893.3245}, 'oxonium': {'ox_HexNAc': 204.0867, 'ox_Hex': 163.0601}}",NNHHH,"[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]"
```


## Encoding Definitions

The following is found in the encoding_definition.txt file, produced using the `-d` flag.

```CSV
Type,Position,Feature,Index,Size
Peptide,1-50,A,0,20
Peptide,1-50,C,1,20
Peptide,1-50,D,2,20
Peptide,1-50,E,3,20
Peptide,1-50,F,4,20
Peptide,1-50,G,5,20
Peptide,1-50,H,6,20
Peptide,1-50,I,7,20
Peptide,1-50,K,8,20
Peptide,1-50,L,9,20
Peptide,1-50,M,10,20
Peptide,1-50,N,11,20
Peptide,1-50,P,12,20
Peptide,1-50,Q,13,20
Peptide,1-50,R,14,20
Peptide,1-50,S,15,20
Peptide,1-50,T,16,20
Peptide,1-50,V,17,20
Peptide,1-50,W,18,20
Peptide,1-50,Y,19,20
Glycan,1-30,N,0,4
Glycan,1-30,H,1,4
Glycan,1-30,F,2,4
Glycan,1-30,A,3,4
Charge,1,Charge_1,0,10
Charge,1,Charge_2,1,10
Charge,1,Charge_3,2,10
Charge,1,Charge_4,3,10
Charge,1,Charge_5,4,10
Charge,1,Charge_6,5,10
Charge,1,Charge_7,6,10
Charge,1,Charge_8,7,10
Charge,1,Charge_9,8,10
Charge,1,Charge_10,9,10
```

## License

This script is released under the MIT License. 

## Acknowledgments

- The Hitchhiker’s Guide to Glycoproteomics.

Oliveira, Tiago, Morten Thaysen-Andersen, Nicolle Packer, and Daniel Kolarich. “The Hitchhiker’s Guide to Glycoproteomics.” Biochemical Society Transactions 49 (July 20, 2021). https://doi.org/10.1042/BST20200879.

- In Silico Platform for Prediction of N-, O- and C-Glycosites in Eukaryotic Protein Sequences.

Chauhan, Jagat Singh, Alka Rao, and Gajendra P. S. Raghava. “In Silico Platform for Prediction of N-, O- and C-Glycosites in Eukaryotic Protein Sequences.” PLoS ONE 8, no. 6 (June 28, 2013): e67008. https://doi.org/10.1371/journal.pone.0067008.

- Large-Scale Identification of N-Linked Intact Glycopeptides in Human Serum Using HILIC Enrichment and Spectral Library Search.

Shu, Qingbo, Mengjie Li, Lian Shu, Zhiwu An, Jifeng Wang, Hao Lv, Ming Yang, et al. “Large-Scale Identification of N-Linked Intact Glycopeptides in Human Serum Using HILIC Enrichment and Spectral Library Search.” Molecular & Cellular Proteomics : MCP 19, no. 4 (April 2020): 672–89. https://doi.org/10.1074/mcp.RA119.001791.

- Assessing the Hydrophobicity of Glycopeptides Using Reversed-Phase Liquid Chromatography and Tandem Mass Spectrometry.

Wang, Junyao, Aiying Yu, Byeong Gwan Cho, and Yehia Mechref. “Assessing the Hydrophobicity of Glycopeptides Using Reversed-Phase Liquid Chromatography and Tandem Mass Spectrometry.” Journal of Chromatography. A 1706 (September 13, 2023): 464237. https://doi.org/10.1016/j.chroma.2023.464237.

- Molecular Basis of C-Mannosylation – a Structural Perspective.

Crine, Samuel L., and K. Ravi Acharya. “Molecular Basis of C-Mannosylation – a Structural Perspective.” The FEBS Journal 289, no. 24 (2022): 7670–87. https://doi.org/10.1111/febs.16265.

- BioPython for handling FASTA files.

Cock, P. J., Antao, T., Chang, J. T., Chapman, B. A., Cox, C. J., Dalke, A., … others. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics, 25(11), 1422–1423.

- Pyteomics for accurate peptide mass calculations.

Goloborodko, A.A.; Levitsky, L.I.; Ivanov, M.V.; and Gorshkov, M.V. (2013) “Pyteomics - a Python Framework for Exploratory Data Analysis and Rapid Software Prototyping in Proteomics”, Journal of The American Society for Mass Spectrometry, 24(2), 301–304. DOI: 10.1007/s13361-012-0516-6

Levitsky, L.I.; Klein, J.; Ivanov, M.V.; and Gorshkov, M.V. (2018) “Pyteomics 4.0: five years of development of a Python proteomics framework”, Journal of Proteome Research. DOI: 10.1021/acs.jproteome.8b00717

- Multicenter Longitudinal Quality Assessment of MS-Based Proteomics in Plasma and Serum.

Kardell, Oliver, Thomas Gronauer, Christine von Toerne, Juliane Merl-Pham, Ann-Christine König, Teresa K. Barth, Julia Mergner, et al. “Multicenter Longitudinal Quality Assessment of MS-Based Proteomics in Plasma and Serum.” Journal of Proteome Research, February 7, 2025. https://doi.org/10.1021/acs.jproteome.4c00644.

- GlyGen: Computational and Informatics Resources for Glycoscience.

York WS, Mazumder R, Ranzinger R, Edwards N, Kahsay R, Aoki-Kinoshita KF, Campbell MP, Cummings RD, Feizi T, Martin M, Natale DA, Packer NH, Woods RJ, Agarwal G, Arpinar S, Bhat S, Blake J, Castro LJG, Fochtman B, Gildersleeve J, Goldman R, Holmes X, Jain V, Kulkarni S, Mahadik R, Mehta A, Mousavi R, Nakarakommula S, Navelkar R, Pattabiraman N, Pierce MJ, Ross K, Vasudev P, Vora J, Williamson T, Zhang W. GlyGen: Computational and Informatics Resources for Glycoscience. Glycobiology. 2020 Jan 28;30(2):72-73. doi: 10.1093/glycob/cwz080. PMID: 31616925; PMCID: PMC7335483.

- Glycosylation of Viral Proteins: Implication in Virus–Host Interaction and Virulence.

Feng, Tingting, Jinyu Zhang, Zhiqian Chen, Wen Pan, Zhengrong Chen, Yongdong Yan, and Jianfeng Dai. “Glycosylation of Viral Proteins: Implication in Virus–Host Interaction and Virulence.” Virulence 13, no. 1 (n.d.): 670–83. https://doi.org/10.1080/21505594.2022.2060464.

- Role of Protein Glycosylation in Interactions of Medically Relevant Fungi with the Host.

Gómez-Gaviria, Manuela, Ana P. Vargas-Macías, Laura C. García-Carnero, Iván Martínez-Duncker, and Héctor M. Mora-Montes. “Role of Protein Glycosylation in Interactions of Medically Relevant Fungi with the Host.” Journal of Fungi 7, no. 10 (October 18, 2021): 875. https://doi.org/10.3390/jof7100875.

-  XMAn v2—a Database of Homo Sapiens Mutated Peptides.

Flores, Marcela Aguilera, and Iulia M Lazar. “XMAn v2—a Database of Homo Sapiens Mutated Peptides.” Bioinformatics 36, no. 4 (February 15, 2020): 1311–13. https://doi.org/10.1093/bioinformatics/btz693.

# Appendix

Additional information, logging runs, and references.

## Log File

The log file provides detailed information about the processing steps, including the number of peptides found after cleavage and the identified N-glycopeptides.

```bash
python glycopeptide_finder_cmd.py -i test_proteomes/human_uniprotkb_proteome_UP000005640_AND_revi_2025_01_17.fasta -p trypsin -c 0 -l log.txt
```

The script generates a log file that records the processing details of each protein sequence. Logging to a text file can be activated with the `-l log.txt` flag. Below are some example log entries:

```log
2025-02-04 00:57:21,942 - INFO - Processing sp|A0A087X1C5|CP2D7_HUMAN with 515 amino acids.
2025-02-04 00:57:21,942 - INFO - Found 50 peptides after trypsin cleavage. The peptides were: ['MGLEALVPLAMIVAIFLLLVDLMHR', 'HQR', 'WAAR', 'YPPGPLPLPGLGNLLHVDFQNTPYCFDQLR', 'R', 'R', 'FGDVFSLQLAWTPVVVLNGLAAVR', 'EAMVTR', 'GEDTADRPPAPIYQVLGFGPR', 'SQGVILSR', 'YGPAWR', 'EQR', 'R', 'FSVSTLR', 'NLGLGK', 'K', 'SLEQWVTEEAACLCAAFADQAGRPFRPNGLLDK', 'AVSNVIASLTCGR', 'R', 'FEYDDPR', 'FLR', 'LLDLAQEGLK', 'EESGFLR', 'EVLNAVPVLPHIPALAGK', 'VLR', 'FQK', 'AFLTQLDELLTEHR', 'MTWDPAQPPR', 'DLTEAFLAK', 'K', 'EK', 'AK', 'GSPESSFNDENLR', 'IVVGNLFLAGMVTTSTTLAWGLLLMILHLDVQR', 'GR', 'R', 'VSPGCPIVGTHVCPVR', 'VQQEIDDVIGQVR', 'RPEMGDQAHMPCTTAVIHEVQHFGDIVPLGVTHMTSR', 'DIEVQGFR', 'IPK', 'GTTLITNLSSVLK', 'DEAVWK', 'KPFR', 'FHPEHFLDAQGHFVKPEAFLPFSAGR', 'R', 'ACLGEPLAR', 'MELFLFFTSLLQHFSFSVAAGQPRPSHSR', 'VVSFLVTPSPYELCAVPR', '']
2025-02-04 00:57:21,943 - INFO - Found 1 N-glycopeptides. The glycopeptides were: [('GTTLITNLSSVLK', 416)]
2025-02-04 00:57:21,943 - INFO - Processing sp|A0A0B4J2F0|PIOS1_HUMAN with 54 amino acids.
2025-02-04 00:57:21,945 - INFO - Found 9 peptides after trypsin cleavage. The peptides were: ['MFR', 'R', 'LTFAQLLFATVLGIAGGVYIFQPVFEQYAK', 'DQK', 'ELK', 'EK', 'MQLVQESEEK', 'K', 'S']
2025-02-04 00:57:21,945 - INFO - Found 0 N-glycopeptides. The glycopeptides were: []
2025-02-04 00:57:21,945 - INFO - Processing sp|A0A0C5B5G6|MOTSC_HUMAN with 16 amino acids.
2025-02-04 00:57:21,945 - INFO - Found 5 peptides after trypsin cleavage. The peptides were: ['MR', 'WQEMGYIFYPR', 'K', 'LR', '']
2025-02-04 00:57:21,945 - INFO - Found 0 N-glycopeptides. The glycopeptides were: []
2025-02-04 00:57:21,945 - INFO - Processing sp|A0A0K2S4Q6|CD3CH_HUMAN with 201 amino acids.
2025-02-04 00:57:21,945 - INFO - Found 13 peptides after trypsin cleavage. The peptides were: ['MTQR', 'AGAAMLPSALLLLCVPGCLTVSGPSTVMGAVGESLSVQCR', 'YEEK', 'YK', 'TFNK', 'YWCR', 'QPCLPIWHEMVETGGSEGVVR', 'SDQVIITDHPGDLTFTVTLENLTADDAGK', 'YR', 'CGIATILQEDGLSGFLPDPFFQVQVLVSSASSTENSVK', 'TPASPTRPSQCQGSLPSSTCFLLLPLLK', 'VPLLLSILGAILWVNRPWR', 'TPWTES']
2025-02-04 00:57:21,945 - INFO - Found 1 N-glycopeptides. The glycopeptides were: [('SDQVIITDHPGDLTFTVTLENLTADDAGK', 100)]
```

### Notes

1. The script assumes well-formatted FASTA input files.
2. Only N-linked glycosylation sequons are detected (no O-linked or other modifications).
3. FASTA protein files contain new lines and or return carrages. When returning to the FASTA, remember this when searching for peptide in original sequence. 

List all common names used in test_proteome folder.

```
for file in ./test_proteomes/*; do
  filename=$(basename "$file")
  part_before_underscore="${filename%%_*}"
  echo "$part_before_underscore"
done
```

## Test Proteomes

Test proteome files from UniProt are available in the `test_proteomes` folder. Below is a list of species gathered. Only Swiss-Prot reviewed proteins were downloaded, and not every sequence available for a species is included. 

I used these test proteomes to generate a zoo of glycopeptides under constrained conditions to fit into a GitHub repo. To build full zoo, remove constraints in batch processing script.

| Common Name   | Scientific Name                                                                 | Taxon ID                |
|---------------|---------------------------------------------------------------------------------|-------------------------|
| Alpaca | Vicugna pacos | 30538 |
| Amoeba | Naegleria gruberi | 5762 |
| Anemone | Nematostella vectensis | 45351 |
| Ant | Camponotus floridanus | 104421 |
| Apple         | Malus domestica                                                                 | 3750                    |
| Arabidopsis   | Arabidopsis thaliana                                                            | 3702                    |
| Aspergillus fumigata | Aspergillus fumigata (strain ATCC MYA-4609 / CBS 101355 / FGSC A1100 / Af293) | 330879 |
| Aspergillus nidulans | Emericella nidulans (strain FGSC A4 / ATCC 38163 / CBS 112.46 / NRRL 194 / M139) | 227321 |
| Avocado | Persea americana | 3435 |
| Banana | Musa acuminata | 4641 | 
| Barley | Hordeum vulgare | 4513 |
| Bat           | Myotis lucifugus                                                                | 59463                   |
| Black Cherry | Prunus serotina | 23207 |
| Black Truffle | Tuber melanosporum (strain Mel28) | 656061 | 
| Blood Fluke | Schistosoma mansoni | 6183 |
| Brine Shrimp | Artemia franciscana | 6661 |
| Brown Alga | Ectocarpus siliculosus | 2880 |
| Bushbaby | Otolemur garnettii | 30611 |
| Camel | Camelus bactrianus | 9837 |
| Candida albicans (Yeast, human pathogen) | Candida albicans (strain SC5314 / ATCC MYA-2876) | 237561 |
| Cat           | Felis catus                                                                     | 9685                    |
| C. elegans    | Caenorhabditis elegans                                                          | 6239                    |
| Chameleon | Anolis carolinensis | 28377 |
| Charcoal Rot | Macrophomina phaseolina (strain MS6) | 1126212 |
| Chicken       | Gallus gallus                                                                   | 9031                    |
| Chimpanzee    | Pan troglodytes                                                                 | 9598                    |
| Chinchilla | Chinchilla lanigera | 34839 |
| C. jejuni     | Campylobacter jejuni                                                            | 1951                    |
| Coffee | Coffea arabica | 13443 |
| Cow           | Bos taurus                                                                      | 9913                    |
| Crocodile | Crocodylus porosus | 8502 |
| Crytococcus | Cryptococcus neoformans var. neoformans serotype D (strain JEC21 / ATCC MYA-565) | 214684|
| Cytomegalovirus | Human cytomegalovirus (strain Merlin) | 295027 |
| Corn Smut | Mycosarcoma maydis | 5270 |
| Date Palm | Phoenix dactylifera | 42345 |
| Debaryomyces hansenii (yeast) | Debaryomyces hansenii (strain ATCC 36239 / CBS 767 / BCRC 21394 / JCM 1990 / NBRC 0083 / IGC 2968) | 284592 |
| Deer Tick | Ixodes scapularis | 6945 | 
| Diatom | Thalassiosira pseudonana | 35128 |
| Dictyostelium | Dictyostelium discoideum                                                        | 44689                   |
| Dog           | Canis lupus familiaris                                                          | 9615                    |
| Donkey        | Equus asinus                                                                    | 9796                    |
| Duck          | Cairina moschata                                                                | 8855                    |
| Dugbe Virus | Dugbe virus (isolate ArD44313) | 766194 |
| Ebola | Zaire ebolavirus (strain Mayinga-76) | 128952 | 
| Elephant      | Loxodonta africana (African Elephant)                                           | 9785                    |
| Fall Armyworm | Spodoptera frugiperda (Fall Armyworm) | 7108 |
| Ferret | Mustela putorius furo | 9669 |
| Fission Yeast | Schizosaccharomyces japonicus (strain yFS275 / FY16936) | 402676 |
| Frog | Xenopus laevis | 8355 |
| Fruit Fly     | Drosophila melanogaster                                                         | 7227                    |
| Goat          | Capra hircus                                                                    | 9925                    |
| Gorilla | Gorilla gorilla gorilla | 9595 |
| Grape | Vitis vinifera | 29760 |
| Green Alga | Chlamydomonas reinhardtii | 3055 |
| Guinea Pig    | Cavia porcellus                                                                 | 10141                   |
| Hamster | Mesocricetus auratus | 10036 | 
| Hemp          | Cannabis sativa                                                                 | 4565                    |
| HHV-1 | Human herpesvirus 1 (strain 17) | 10299 | 
| HIV-1 | Human immunodeficiency virus type 1 group N (isolate YBF30) | 388818 |
| HIV-2 | Human immunodeficiency virus type 2 subtype A (isolate BEN) | 11714 | 
| Honeybee      | Apis mellifera                                                                  | 7460                    |
| Horse         | Equus caballus                                                                  | 9796                    |
| HRSV S-2 | Human respiratory syncytial virus A (strain S-2) | 410078 | 
| Human         | Homo sapiens                                                                    | 9606                    |
| Influenza B | Influenza B virus (strain B/Lee/1940) | 518987 |
| Influenza C | Influenza C virus (strain C/Ann Arbor/1/1950) | 11553 |
| JEV | Japanese encephalitis virus (strain M28) | 2555554 |
| Kidney Bean | Phaseolus vulgaris | 3885 |
| Kluyveromyces lactis (lactate processing yeast) | Kluyveromyces lactis (strain ATCC 8585 / CBS 2359 / DSM 70799 / NBRC 1267 / NRRL Y-1140 / WM37) | 284590 |
| LASV | Lassa virus (strain Mouse/Sierra Leone/Josiah/1976) | 11622 |
| LCMV | Lymphocytic choriomeningitis virus (strain Armstrong) | 11624 |
| Lemur | Microcebus murinus | 30608 |
| Macaque (Rhesus monkey) | Macaca mulatta | 9544 |
| Maize | Zea mays | 4577 |
| Measles virus | Measles virus (strain Ichinose-B95a) | 645098 |
| Monkey (cynomolgus, crab-eating) | Macaca fascicularis | 9541 |
| Mosquito (African malaria) | Anopheles gambiae | 7165 |
| Mouse         | Mus musculus                                                                    | 10090                   |
| Naked Mole Rat | Heterocephalus glaber | 10181 |
| Nematode (roundworm) | Caenorhabditis briggsae | 6238 |
| Norovirus | Norovirus (strain Human/NoV/United States/Norwalk/1968/GI) | 524364 |
| Octopus | Octopus vulgaris | 6645 |
| Olive | Olea europaea | 4146 |
| Opossum | Monodelphis domestica | 13616 | 
| Orange | Citrus sinensis | 2711 |
| Orangutan     | Pongo abelii                                                                    | 9601                    |
| Oyster | Magallana gigas | 29159 |
| Paramecium | Paramecium tetraurelia | 5888 |
| Peach | Prunus persica | 3760 | 
| Penicillium | Penicillium rubens (strain ATCC 28089 / DSM 1075 / NRRL 1951 / Wisconsin 54-1255) | 500485 | 
| Pig (Domestic)           | Sus scrofa domesticus   | 9823                    |
| Platypus | Ornithorhynchus anatinus | 9258 | 
| Poplar Leaf Rust Fungus | Melampsora larici-populina (strain 98AG31 / pathotype 3-4-7) | 747676 |
| Potato | Solanum tuberosum | 4113 |
| Psilocybe mushroom | Psilocybe cubensis| 181762 |
| Pufferfish | Takifugu rubripes | 31033 | 
| Rabbit | Oryctolagus cuniculus | 9986 |
| Rat           | Rattus norvegicus                                                               | 10116                   |
| Red Alga | Cyanidioschyzon merolae (strain NIES-3377 / 10D) | 280699 |
| Rice          | Oryza sativa subsp. japonica                                                    | 39947                   |
| Rice Blast Fungus | Pyricularia oryzae (strain 70-15 / ATCC MYA-4617 / FGSC 8958) | 242507 |
| Rice Fish (Japanese) | Oryzias latipes | 8090 | 
| RVA | Rotavirus A (isolate RVA/Monkey/South Africa/SA11-H96/1958/G3P5B[2]) | 450149 |
| RVB | Rotavirus B (isolate RVB/Human/China/ADRV/1982) | 10942 |
| RVC | Rotavirus C (isolate RVC/Human/United Kingdom/Bristol/1989) | 31567 |
| SARS-CoV      | SARS-CoV (Severe Acute Respiratory Syndrome Coronavirus)                        | 694009                  |
| SFTSV | SFTS phlebovirus (isolate SFTSV/Human/China/HB29/2010) | 992212 |
| Shark | Callorhinchus milii | 7868 | 
| Sheep         | Ovis aries                                                                      | 9940                    |
| Silk Moth | Bombyx mori | 7091 |
| Silveira (Coccidioides Silveira strain) | Coccidioides posadasii (strain RMSCC 757 / Silveira) | 443226 | 
| Snake (Brown Eastern) | Pseudonaja textilis | 8673 | 
| Softshell Turtle | Pelodiscus sinensis | 13735 | 
| Spike Moss (lycophyte) | Selaginella moellendorffii | 88036 | 
| Sponge | Amphimedon queenslandica | 400682 |
| Sorghum       | Sorghum bicolor                                                                 | 4558                    |
| Squirrel      | Ictidomys tridecemlineatus                                                      | 43179                   |
| Starfish | Patiria pectinifera | 7594 |
| Strawberry | Fragaria ananassa | 3747 |
| Sugarcane | Saccharum officinarum | 4547 |
| Sunflower | Helianthus annuus | 4232 |
| Sycamore | Platanus occidentalis | 4403 |
| Tea plant | Camellia sinensis | 4442 |
| Tobacco | Nicotiana tabacum  | 4097 |
| Tilapia | Oreochromis niloticus | 8128 | 
| Tomato | Solanum lycopersicum | 4081 |
| Trout (Brown) | Oreochromis niloticus | 8128 | 
| Turkey | Meleagris gallopavo | 9103 | 
| Urchin | Strongylocentrotus purpuratus | 7668 |
| VZV | Varicella-zoster virus (strain Dumas) | 10338 |
| Wasp (parasitoid) | Nasonia vitripennis | 7425 |
| Watermelon | Citrullus lanatus | 3654 |
| Wheat         | Triticum aestivum                                                               | 4565                    |
| Whisk fern | Psilotum nudum | 3240 |
| Wild Rice (North America) | Oryza nivara | 4536 | 
| WNV | West Nile virus | 11082 | 
| XMAn v2 Missense | Homo sapians - Unknown Mutation Analysis (Human missense peptide library) Download at: https://github.com/lazarlab/XMAn-v2 | 9606 |
| XMAn v2 Nonsense| Homo sapians - Unknown Mutation Analysis (Human nonsense peptide library) Download at: https://github.com/lazarlab/XMAn-v2 | 9606 |
| Yak | Bos mutus grunniens | 30521 |
| Yeast (Budding, Baker's)        | Saccharomyces cerevisiae (strain ATCC 204508 / S288c)                           | 559292                  |
| Yeast (Fission) | Schizosaccharomyces pombe (strain 972 / ATCC 24843) | 284812 |
| Zebra Finch | Taeniopygia guttata | 59729 | 
| Zebrafish     | Danio rerio                                                                     | 7955                    |
| Zebu | Bos indicus | 9915 | 
| Zika | Zika virus | 64320 | 

## Glycan Mass Library

Glycan mass libraries were gathered from GlyGen. The follow data was processed and used in this tool. Plug and play glycans to meet needs.

It is advisable to create a targeted glycan library along with a list of peptides to compute an intact glycan library. The size of the library can grow rapidly, so it is important to manage it effectively.

| File Name          | Glycan Count |
|--------------------|--------------|
| glycan_database     | 44686        |
| glycan_type_n_linked_byonic.csv | 369    |
| test_glycan_library.csv | 3       |
