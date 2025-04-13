"""
glycopeptide_sequence_finder_cmd.py

This script processes a FASTA file to identify glycopeptides based on protease cleavage rules and glycosylation sequons, 
predicts their masses, calculates their hydrophobicity, and calculates their pI. The results are written to a CSV file 
in the digest_glycopeptide_library directory.

Usage:
    python glycopeptide_sequence_finder_cmd.py -i <input_fasta_file> -o <output_csv_file> -p <protease> -g <glycosylation_type> -c <missed_cleavages> -l <log_file> -v -y <glycan_file> -z <max_charge> --ion-series

Author:
    Richard Shipman -- 2025
"""
import argparse
import csv
from itertools import permutations
import re
from Bio import SeqIO
import os
import logging
import pandas as pd

# Constants

# Define protease cleavage rules
proteases = {
    "trypsin": ("[KR]", "P"),  # Cleaves after K or R unless followed by P
    "chymotrypsin": ("[FLWY]", "P"),  # Cleaves after F, L, W, or Y unless followed by P
    "glu-c": ("[DE]", "P"),  # Cleaves after E
    "lys-c": ("K", "P"),  # Cleaves after K unless followed by P
    "arg-c": ("R", "P"),  # Cleaves after R
    "pepsin": ("[FWY]", None),  # Cleaves after F, W, or Y
    "asp-n": ("D", None),  # Cleaves **before** Asp (D)
    "proteinase-k": ("[AFILVWY]", None),  # Cleaves after A, F, I, L, V, W, Y
    "operator" : ("[ST]", None), # Cleaves after glycosylated S or T # set missed cleavage to 10 # for O glycopeptides
    "thermolysin": ("[ALIVFM]", "[DE]") # Cleaves after A, L, I, V, F, M unless followed by D or E
}

# Define glycosylation sequon rules
glycosylation = {
    "N": ("N[^P][ST]"),  # N-glycosylation sequon - N-X-S/T (X is any amino acid except P).
    "O": ("[ST]"),  # O-glycosylation sequon - S/T (experimental! Creates large number of O-glycopeptides.)
    "C": ("W..[WCF]"),  # C-glycosylation sequon - W-X-X-W, W-X-X-C, or W-X-X-F (X is any amino acid).

    # Additional motifs
    "C-Mannose-1": ("W..W..W"),  # W-X-X-W-X-X-W
    "C-Mannose-2": ("W..W..W..C"),  # W-X-X-W-X-X-W-X-X-C
    "C-Mannose-3": ("W..[WC]"),  # W-X-X-W or W-X-X-C
    "N-Glycan": ("N[^P][STC]"),  # Extended N-X-S/T/C
    "O-Fucose-1": ("C....[ST]C"),  # C-X-X-X-X-S/T-C
    "O-Fucose-2": ("C..[ST]C"),  # C-X-X-S/T-C
    "O-Fucose-3": ("C..[ST]C..G"),  # C-X-X-S/T-C-X-X-G
    "O-GlcNAc": ("C....[ST]G..C"),  # C-X-X-X-X-S/T-G-X-X-C
    "O-Glucose-1": ("C.NT.GS[FY].C"),  # C-X-N-T-X-G-S-(F/Y)-X-C (corrected from previous entry)
    "O-Glucose-2": ("C.S.[PA]C"),  # C-X-S-X-P/A-C

    # mutation glycosylation site motifs
    "type-mutation-c": ("[^N][^P][ST]"),  # N-glycosylation mutation at first position
    "type-mutation-b": ("NP[ST]"),  # N-glycosylation mutation at second pos
    "type-mutation-a": ("N[^P][^S^T]"),  # N-glycosylation mutation at third pos
    "type-extended-mutation-a": ("[^N][^P][STC]"),  # N-glycosylation mutation at first position
    "type-extended-mutation-b": ("NP[STC]"),  # N-glycosylation mutation at second pos
    "type-extended-mutation-c": ("N[^P][^S^T^C]"),  # N-glycosylation mutation at third pos extended to include C

}

# Define amino acid mass values for calculating peptide mass
amino_acid_masses = {
    'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694, 'C': 103.00919,
    'Q': 128.05858, 'E': 129.04259, 'G': 57.02146, 'H': 137.05891, 'I': 113.08406,
    'L': 113.08406, 'K': 128.09496, 'M': 131.04049, 'F': 147.06841, 'P': 97.05276,
    'S': 87.03203, 'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
}

# Kyte-Doolittle hydrophobicity values for amino acids
hydrophobicity_values = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# Define amino acid pKa values for calculating pI (placeholder)
pKa_values = {
    'C': 8.18, 'D': 3.65, 'E': 4.25, 'H': 6.00, 'K': 10.53,
    'R': 12.48, 'Y': 10.07, 'N': 3.22, 'Q': 3.22, 'S': 3.70,
    'T': 3.70, 'W': 10.07
}

# Define default N-glycan mass library as a DataFrame (N-Glycans) Using most common N-glycans as default, HexNAc(2)Hex(8) with no stereochemistry
default_n_glycan_library = pd.DataFrame([
    #{"glytoucan_ac": "G59324HL", "byonic": "HexNAc(2)Hex(12) % 2350.792627", "composition": "HexNAc(2)Hex(12)", "mass": 2350.792627, "shorthand_glycan": "N2H12"}, # N2H12
    #{"glytoucan_ac": "G58087IP", "byonic": "HexNAc(2)Hex(11) % 2188.739804", "composition": "HexNAc(2)Hex(11)", "mass": 2188.739804, "shorthand_glycan": "N2H11"}, # N2H11
    #{"glytoucan_ac": "G83460ZZ", "byonic": "HexNAc(2)Hex(10) % 2026.686980", "composition": "HexNAc(2)Hex(10)", "mass": 2026.686980, "shorthand_glycan": "N2H10"}, # N2H10
    #{"glytoucan_ac": "G80920RR", "byonic": "HexNAc(2)Hex(9) % 1864.634157", "composition": "HexNAc(2)Hex(9)", "mass": 1864.634157, "shorthand_glycan": "N2H9"}, # N2H9
    #{"glytoucan_ac": "G62765YT", "byonic": "HexNAc(2)Hex(8) % 1702.581333", "composition": "HexNAc(2)Hex(8)", "mass": 1702.581333, "shorthand_glycan": "N2H8"}, # N2H8 -- High Mannose
    #{"glytoucan_ac": "G31852PQ", "byonic": "HexNAc(2)Hex(7) % 1540.528510", "composition": "HexNAc(2)Hex(7)", "mass": 1540.528510, "shorthand_glycan": "N2H7"}, # N2H7
    #{"glytoucan_ac": "G41247ZX", "byonic": "HexNAc(2)Hex(6) % 1378.475686", "composition": "HexNAc(2)Hex(6)", "mass": 1378.475686, "shorthand_glycan": "N2H6"}, # N2H6
    {"glytoucan_ac": "G22768VO", "byonic": "HexNAc(2)Hex(3) % 1216.422863", "composition": "HexNAc(2)Hex(3)", "mass": 1216.422863, "shorthand_glycan": "N2H3"}, # N2H3
    #{"glytoucan_ac": "G36670VW", "byonic": "HexNAc(5)Hex(5)dHex(1)NeuAc(2) % 2553.909723", "composition": "HexNAc(5)Hex(5)dHex(1)NeuAc(2)", "mass": 2553.909723, "shorthand_glycan": "N5H5F1S2"}  # N5H5F1S2 -- # Complex - Fucosylation, Sialylated
    #{"glytoucan_ac": "G29068FM", "byonic": "HexNAc(1) % 221.089937305", "composition": "HexNAc(1)", "mass": 221.089937305, "shorthand_glycan": "N1"}, # N1 -- EndoH Treatment HexNAc
    #{"glytoucan_ac": "G04038LG", "byonic": "HexNAc(1)dHex(1) % 367.147846175", "composition": "HexNAc(1)dHex(1)", "mass": 367.147846175, "shorthand_glycan": "N1F1"}, # N1F1 -- EndoH Treatment HexNAc (Fuc-Core)
])

# Define default O-glycan mass library as a DataFrame (O-GalNAc Glycans) Using most common O-glycans as default, HexNac(1) with no stereochemistry
default_o_glycan_library = pd.DataFrame([
    {"glytoucan_ac": "G14843DJ", "byonic": "HexNAc(1) % 221.089937305", "composition": "HexNAc(1)", "mass": 221.089937305, "shorthand_glycan": "N1"}, # N1
])

# Define default C-glycan mass library as a DataFrame (C-Mannosylation) Using most common C-glycans as default, Hex(1) with no stereochemistry
default_c_glycan_library = pd.DataFrame([
    {"glytoucan_ac": "G81399MY", "byonic": "Hex(1) % 180.0633882", "composition": "Hex(1)", "mass": 180.0633882, "shorthand_glycan": "H1"}, # H1
])

# Define monosaccharide mass library with multiple properties (later use)
monosaccharide_library = {
    "Hex": {"mass": 162.0528, "formula": "C6H10O5", "symbol": "H"},
    "HexA": {"mass": 176.0321, "formula": "C6H8O6", "symbol": "Ha"},
    "HexNAc": {"mass": 203.0794, "formula": "C8H13NO5", "symbol": "N"},
    "Fuc": {"mass": 146.0579, "formula": "C6H12O5", "symbol": "F"},
    "dHex": {"mass": 146.0579, "formula": "C6H12O5", "symbol": "dH"},
    "NeuAc": {"mass": 291.0954, "formula": "C11H17NO8", "symbol": "S"},
    "NeuGc": {"mass": 307.0903, "formula": "C11H17NO9", "symbol": "G"},
    "Xyl": {"mass": 150.0423, "formula": "C5H10O5", "symbol": "X"},
    "Pent": {"mass": 132.0423, "formula": "C5H10O5", "symbol": "P"},
    "Sulpho": {"mass": 79.9568, "formula": "SO3", "symbol": "Su"},
    "Phospho": {"mass": 79.9663, "formula": "PO3", "symbol": "Ph"},
    "Methyl": {"mass": 14.0157, "formula": "CH3", "symbol": "Me"},
    "Acetyl": {"mass": 42.0106, "formula": "C2H3O", "symbol": "Ac"},
    "Deoxy": {"mass": -18.0106, "formula": "H2O", "symbol": "d"},
    "Amine": {"mass": 1.0078, "formula": "H", "symbol": "NH2"}
}

# Define hydrophobicity values for glycan library (experimental) 
# (between -1 and 1) in relation to the peptide backbone and rank spread of glycopeptides with the same backdone.
# compute_glycan_hydrophobicity.py script will be used to calculate the hydrophobicity values.
# averaged aggregated from human and mouse datasets
default_glycan_hydrophobicity = {
    "N2H8": -0.0062899,
    "N2H7": -0.005723031,
    "N2H6": -0.004347122,
    "N2H9": -0.001982055,
    "N2H10": -0.001554309,
    "N2H11": -2.53305E-05,
    "N2H12": 5.36847E-06
}

# Add these constants after the existing mass constants
# Immonium ion masses (mass of residue minus CO)
immonium_masses = {
    'A': 44.05003, 'R': 129.11400, 'N': 87.05584, 'D': 88.03986, 'C': 76.02210,
    'E': 102.05550, 'Q': 101.07150, 'G': 30.03438, 'H': 110.07180, 'I': 86.09698,
    'L': 86.09698, 'K': 101.10790, 'M': 104.05340, 'F': 120.08130, 'P': 70.06568,
    'S': 60.04494, 'T': 74.06059, 'W': 159.09220, 'Y': 136.07620, 'V': 72.08133
}

# Functions

def cleave_sequence(sequence, protease, missed_cleavages=0):
    """
    Cleaves a sequence based on protease rules and filters by minimum length.
    Ensures no duplicate peptides are generated when using missed cleavages.
    
    Parameters:
        sequence (str): The protein sequence to cleave
        protease (str): The protease to use for cleavage
        missed_cleavages (int): Number of missed cleavages allowed (default: 0)
    
    Returns:
        list: List of unique peptides meeting the length criteria
    """
    cleavage_pattern, exclusion = proteases[protease.lower()]

    # Handle Asp-N separately because it cleaves **before** D
    if protease.lower() == "asp-n":
        regex = rf"(?={cleavage_pattern})"  # Cleaves before D
    else:
        regex = rf"(?<={cleavage_pattern})(?!{exclusion})"  # Cleaves after other residues

    # Perform cleavage
    fragments = re.split(regex, sequence)

    # Use a set to store unique peptides
    unique_peptides = set()
    
    # Generate peptides including missed cleavages
    for i in range(len(fragments)):
        for j in range(i + 1, min(i + 2 + missed_cleavages, len(fragments) + 1)):
            peptide = "".join(fragments[i:j])
            # Only include peptides with length >= 5
            if len(peptide) >= 5:
                unique_peptides.add(peptide)

    # Convert set back to sorted list for consistent output
    return sorted(list(unique_peptides))

def find_glycopeptides(peptides_df, glycosylation_type):
    """Identifies peptides containing sequons and maps the sites to the full protein sequence."""
    
    # Process pandas Dataframe
    peptides = [pep for sublist in peptides_df["Peptides"].tolist() for pep in sublist]
    full_sequence = peptides_df["Sequence"].tolist()[0]
    glycosylation_type = glycosylation_type
    
    glyco_sequon = re.compile(glycosylation[glycosylation_type])
    glycopeptides = []
    for pep in peptides:
        for match in glyco_sequon.finditer(pep):
            # Map the position in the peptide to the full protein sequence
            start_in_protein = full_sequence.find(pep) + match.start() + 1  # Adjust to 1-based indexing
            glycopeptides.append((pep, start_in_protein))
    return glycopeptides

def calculate_peptide_mass(sequence):
    """Calculates the mass of a peptide using predefined amino acid masses."""
    
    # Common ambiguous residues
    invalid_residues = {"X", "B", "Z", "J", "U", "O"}  
    if any(aa in invalid_residues for aa in sequence):
        return "Unknown"  # Or return unknown if you prefer
    
    # water
    water  = 18.010565

    # Calculate the mass using the amino_acid_masses dictionary
    mass = sum(amino_acid_masses.get(aa, 0) for aa in sequence) + water
    return mass

def predict_hydrophobicity(peptide_sequence):
    """Predicts the hydrophobicity of a peptide sequence using the Kyte-Doolittle scale."""
    
    # safety check for empty peptide strings
    if len(peptide_sequence) == 0:
        return 0.0

    # Calculate the average hydrophobicity of the peptide
    total_hydrophobicity = sum(hydrophobicity_values.get(aa, 0) for aa in peptide_sequence)
    average_hydrophobicity = round(total_hydrophobicity / len(peptide_sequence), 5)

    return average_hydrophobicity

def calculate_pI(peptide_sequence):
    """
    Calculate the isoelectric point (pI) of a peptide sequence.
    
    Parameters:
        peptide (str): Peptide sequence (1-letter amino acid codes)
    
    Returns:
        float: Estimated isoelectric point (pI)
    """
    # peptide sequence to string
    peptide_sequence = str(peptide_sequence)

    # Terminal group pKa values (assuming N-term = 9.6, C-term = 2.3)
    N_term_pKa = 9.6
    C_term_pKa = 2.3

    # Count occurrences of ionizable residues
    residue_counts = {aa: peptide_sequence.count(aa) for aa in pKa_values}
    
    # Function to calculate charge at a given pH
    def calculate_net_charge(pH):
        charge = 0.0

        # N-terminal charge
        charge += 1 / (1 + 10**(pH - N_term_pKa))

        # C-terminal charge
        charge -= 1 / (1 + 10**(C_term_pKa - pH))

        # Side chain charges
        for aa, count in residue_counts.items():
            if count > 0:
                pKa = pKa_values[aa]
                if aa in ['D', 'E', 'Y', 'C']:  # Acidic side chains
                    charge -= count / (1 + 10**(pKa - pH))
                elif aa in ['H', 'K', 'R']:  # Basic side chains
                    charge += count / (1 + 10**(pH - pKa))
        
        return charge

    # Use bisection method to find the pH where net charge is closest to zero
    low, high = 0.0, 14.0
    while high - low > 0.01:  # Precision threshold
        mid = (low + high) / 2
        net_charge = calculate_net_charge(mid)
        if net_charge > 0:
            low = mid
        else:
            high = mid

    return round((low + high) / 2, 2)

def compute_mz(mass, charge):
    """Compute m/z value for a given mass and charge state."""
    proton = 1.007276
    return (mass + (charge * proton)) / charge

# experimental glycopeptide hydrophobicity calculation
def compute_hf_experimental(peptide_hydrophobicity, glycan, hf_weight=10, rt_scale=60):
    """Compute HF_experimental and scaled retention time (rt_HF_experimental) values."""
    glycan_hydrophobicity = default_glycan_hydrophobicity.get(glycan, None)
    
    # Compute HF_experimental and scaled retention time (rt_HF_experimental) values
    if glycan_hydrophobicity is not None:
        hf_experimental = peptide_hydrophobicity + glycan_hydrophobicity * hf_weight
        rt_hf_experimental = (glycan_hydrophobicity * hf_weight + 1) * (rt_scale / 2)
    else:
        hf_experimental = ""
        rt_hf_experimental = ""
    
    return round(hf_experimental, 5), rt_hf_experimental

def process_glycopeptides(peptide_file, glycans, max_charge):
    """
    Generate glycopeptide variants by combining peptides and glycans, and compute their mass-to-charge (m/z) values.
    Parameters:
        peptide_file (str): Path to the CSV file containing peptide data. The CSV file should include columns such as
                            'ProteinID', 'Site', 'Peptide', 'PredictedMass', 'Length', 'Sequon', 'Hydrophobicity', and 'pI'.
        glycans (pandas.DataFrame): DataFrame containing glycan data. This DataFrame should include columns such as
                                    'glytoucan_ac', 'composition', and 'mass'.
        max_charge (int): The maximum charge state for which m/z values will be computed. m/z values are calculated for charge
                          states ranging from 2 up to and including max_charge.
    Returns:
        pandas.DataFrame: A DataFrame where each row represents a glycopeptide, including the following columns:
    Notes:
        - The peptide CSV file is loaded with low_memory=False to improve type inference.
        - Non-numeric values in the 'PredictedMass' column are coerced to NaN and subsequently dropped.
        - It is assumed that the compute_mz() function is defined in the same scope where process_glycopeptides() is used.
    """
    
    # Load peptide and glycan data
    peptides = pd.read_csv(peptide_file, low_memory=False)
    
    # Convert relevant columns to float
    peptides = peptides[pd.to_numeric(peptides['PredictedMass'], errors='coerce').notnull()]
    peptides['PredictedMass'] = peptides['PredictedMass'].astype(float)
    glycans['mass'] = glycans['mass'].astype(float)
    
    # Initialize list to store results
    results = []
    
    # Generate glycopeptides, compute m/z values, compute HF_experimental
    for _, pep in peptides.iterrows():
        for _, gly in glycans.iterrows():

            # Compute glycopeptide mass
            glycopeptide_mass = pep['PredictedMass'] + gly['mass']

            # Compute m/z values for charge states from 2 to max_charge
            mz_values = {f'z{z}': compute_mz(glycopeptide_mass, z) for z in range(2, max_charge + 1)}

            # Create a dictionary with the results
            result = {
                'ProteinID': pep['ProteinID'],
                'Site': pep['Site'],
                'GlyToucan_AC': gly['glytoucan_ac'],
                'Composition': gly['composition'],
                'ShorthandGlycan': gly['shorthand_glycan'],
                'Peptide': pep['Peptide'],
                'Start': pep['Start'],
                'End': pep['End'],
                'Length': pep['Length'],
                'Sequon': pep['Sequon'],
                'GlycopeptideMass': glycopeptide_mass,
                'PeptideMass': pep['PredictedMass'],
                'GlycanMass': gly['mass'],
                'Hydrophobicity': pep['Hydrophobicity'],
                'pI': pep['pI'],
                **mz_values, # z charge states values
            }

            results.append(result)
    
    return pd.DataFrame(results)

def setup_logging(log_file):
    """Sets up logging to a file."""
    logging.basicConfig(filename=log_file, level=logging.INFO, 
                        format='%(asctime)s - %(levelname)s - %(message)s')

def process_fasta(file, protease, missed_cleavages, glycosylation_type):
    """Processes the input FASTA file and extracts glycopeptides."""
    
    # Initialize list to store results
    results = []

    # Regular expression to capture OS and OX from the header
    os_ox_pattern = r"OS=([^\s]+(?: [^\s]+)*)\s+OX=(\d+)\s+GN=([^\s]+)\s+PE=(\d+)\s+SV=(\d+)"

    # Process each record in the FASTA file
    for record in SeqIO.parse(file, "fasta"):
        protein_id = record.id
        sequence = str(record.seq)

        # Get the full description line
        header = record.description  
        
        # Use regex to find OS, OX, GN, PE, and SV
        match = re.search(os_ox_pattern, header)
        if match:
            os_value = match.group(1)  # Organism Species
            ox_value = match.group(2)  # Organism Taxonomy ID
            gn_value = match.group(3)  # Gene Name
            pe_value = match.group(4)  # Protein Evidence
            sv_value = match.group(5)  # Sequence Version
        else:
            # blank values if not found
            os_value = ""
            ox_value = ""
            gn_value = ""
            pe_value = ""
            sv_value = ""

        # Log the digestion processing of the protein with protease and missed cleavages
        logging.info(f"Processing {protein_id} with {len(sequence)} amino acids.")        
        peptides = cleave_sequence(sequence, protease, missed_cleavages)
        logging.info(f"Found {len(peptides)} peptides after {protease} cleavage. The peptides were: {peptides}")
        
        # Write protein, peptide and data above list of string results to a pandas DataFrame
        results.append({
            "ProteinID": protein_id,
            "Peptides": peptides,
            "Sequence": sequence,
            "Protease": protease,
            "MissedCleavages": missed_cleavages,
            "GlycosylationType": glycosylation_type,
            "Species": os_value,
            "TaxonID": ox_value,
            "GeneName": gn_value,
            "ProteinEvidence": pe_value,
            "SequenceVersion": sv_value
        })

    return pd.DataFrame(results)

def calculate_immonium_ions(peptide):
    """Calculate immonium ions for a peptide sequence."""
    # Get unique amino acids in the peptide
    unique_aas = set(peptide)
    
    # Calculate immonium ions for each unique amino acid
    immonium_ions = {}
    for aa in unique_aas:
        if aa in immonium_masses:
            immonium_ions[f'im_{aa}'] = round(immonium_masses[aa], 4)
    
    return dict(sorted(immonium_ions.items(), key=lambda x: x[1]))

def calculate_internal_fragments(peptide, charge=1):
    """Calculate internal fragment ions for a peptide sequence."""
    proton = 1.007276
    water = 18.010565
    
    internal_ions = {}
    
    # Generate all possible internal fragments of length 2 or more
    for i in range(len(peptide)-1):
        for j in range(i+2, len(peptide)+1):
            fragment = peptide[i:j]
            if len(fragment) >= 2:  # Only consider fragments of length 2 or more
                # Calculate mass of internal fragment
                mass = sum(amino_acid_masses[aa] for aa in fragment)
                # Add proton for charge
                mz = (mass + proton) / charge
                internal_ions[f'int_{i+1}_{j+1}'] = round(mz, 4)
    
    return dict(sorted(internal_ions.items(), key=lambda x: x[1]))

def calculate_peptide_ions(peptide, charge=1):
    """
    Calculate theoretical m/z values for b, y, c, z ions, immonium ions, and internal fragments for a peptide.
    """
    # Get the existing ion calculations
    ions = {
        'b': [],
        'y': [],
        'c': [],
        'z': []
    }
    
    # Constants for ion calculations
    proton = 1.007276
    water = 18.010565
    NH3 = 17.0265

    # Calculate b ions (N-terminal fragments)
    cumulative = 0.0
    for i in range(len(peptide) - 1):
        cumulative += amino_acid_masses[peptide[i]]
        ions['b'].append(round((cumulative + proton) / charge, 4))
    ions['b'].sort()

    # Calculate y ions (C-terminal fragments)
    cumulative = 0.0
    for i in range(len(peptide) - 1, 0, -1):
        cumulative += amino_acid_masses[peptide[i]]
        ions['y'].append(round((cumulative + water + proton) / charge, 4))
    ions['y'].sort()

    # Calculate c ions (N-terminal fragments + NH3)
    cumulative = 0.0
    for i in range(len(peptide) - 1):
        cumulative += amino_acid_masses[peptide[i]]
        ions['c'].append(round((cumulative + NH3 + proton) / charge, 4))
    ions['c'].sort()

    # Calculate z ions (C-terminal fragments - NH3)
    cumulative = 0.0
    for i in range(len(peptide) - 1, 0, -1):
        cumulative += amino_acid_masses[peptide[i]]
        ions['z'].append(round((cumulative + proton - NH3) / charge, 4))
    ions['z'].sort()

    # Add immonium ions
    ions['immonium'] = calculate_immonium_ions(peptide)
    
    # Add internal fragments
    ions['internal'] = calculate_internal_fragments(peptide, charge)

    return ions

def calculate_n_glycopeptide_ions(peptide, glycan_composition, glycan_frag_order=None, charge=1):
    """
    Calculate theoretical m/z values for b, y, c, z, Y, B, oxonium ions, immonium ions, internal ions,
    and Y ions with peptide losses for an N-glycopeptide.
    All ion series are sorted from smallest to largest m/z.

    Parameters:
      peptide (str): The peptide sequence (e.g., "NTSK").
      glycan_composition (str): A string like "HexNAc(5)Hex(5)dHex(1)NeuAc(2)".
      glycan_frag_order (list, optional): A list specifying the sugar loss order.
      charge (int): The charge state (default is 1).

    Returns:
      dict: A dictionary with keys for ion types and their m/z values.
    """
    # Get peptide backbone ions first
    peptide_ions = calculate_peptide_ions(peptide, charge)
    
    # Constants
    proton = 1.007276
    water = 18.010565

    # Find the glycosite position (N-X-S/T motif)
    glycosite_match = re.search(r'N[^P][ST]', peptide)
    if not glycosite_match:
        return peptide_ions  # Return only peptide ions if no glycosite found
    
    glycosite_start = glycosite_match.start()
    glycosite_end = glycosite_match.end()

    # Compute the neutral mass of the peptide (including water)
    peptide_mass = sum(amino_acid_masses[aa] for aa in peptide) + water

    # Parse glycan composition into a dictionary
    glycan_dict = {}
    for part in glycan_composition.split(')'):
        if part:
            sugar, count = part.split('(')
            glycan_dict[sugar] = int(count)

    # --- Calculate immonium ions ---
    immonium_ions = {}
    for aa in peptide:
        if aa in immonium_masses:
            immonium_ions[f'im_{aa}'] = round(immonium_masses[aa] + proton, 4)

    # --- Calculate internal ions ---
    internal_ions = {}
    for i in range(len(peptide) - 2):
        for j in range(i + 2, len(peptide)):
            # Skip if internal fragment would break the glycosite
            if i < glycosite_start and j > glycosite_end:
                continue
            internal_seq = peptide[i:j+1]
            internal_mass = sum(amino_acid_masses[aa] for aa in internal_seq)
            internal_ions[f'internal_{i+1}_{j+1}'] = round((internal_mass + proton) / charge, 4)

    # --- Calculate Y ions (glycan-attached peptide fragments) ---
    Y_ions = {}
    current_mass = peptide_mass
    Y_ions['Y0'] = round((current_mass + proton) / charge, 4)

    # Process glycan additions in order
    remaining_sugars = []
    for sugar, count in glycan_dict.items():
        remaining_sugars.extend([sugar] * count)

    # Sort sugars by mass for consistent ordering
    remaining_sugars.sort(key=lambda x: monosaccharide_library[x]['mass'])

    # Add sugars in order
    for i, sugar in enumerate(remaining_sugars, start=1):
        current_mass += monosaccharide_library[sugar]['mass']
        Y_ions[f'Y{i}'] = round((current_mass + proton) / charge, 4)

    # --- Calculate Y ions with peptide backbone losses ---
    Y_peptide_loss = {}
    
    # Calculate b-ion masses (without charge) - only for positions that retain glycosite
    b_masses = []
    cumulative = 0.0
    for i in range(len(peptide) - 1):
        cumulative += amino_acid_masses[peptide[i]]
        # Only include b-ions that would leave the glycosite intact
        if i < glycosite_start:
            b_masses.append((i+1, cumulative))

    # Calculate y-ion masses (without charge) - only for positions that retain glycosite
    y_masses = []
    cumulative = 0.0
    for i in range(len(peptide) - 1, 0, -1):
        cumulative += amino_acid_masses[peptide[i]]
        # Only include y-ions that would leave the glycosite intact
        if i > glycosite_end:
            y_masses.append((len(peptide) - i, cumulative + water))

    # For each Y ion, subtract each valid b-ion and y-ion mass
    for y_num, y_mass in Y_ions.items():
        y_neutral = (y_mass * charge) - proton
        
        # Subtract b-ions (N-terminal losses)
        for b_idx, b_mass in b_masses:
            fragment_mass = y_neutral - b_mass
            Y_peptide_loss[f'{y_num}-b{b_idx}'] = round((fragment_mass + proton) / charge, 4)
        
        # Subtract y-ions (C-terminal losses)
        for y_idx, y_mass in y_masses:
            fragment_mass = y_neutral - y_mass
            Y_peptide_loss[f'{y_num}-y{y_idx}'] = round((fragment_mass + proton) / charge, 4)

    # --- Calculate B ions (glycan fragment ions) ---
    B_ions = {}
    cumulative = 0.0
    
    # Sort sugars by mass
    sorted_sugars = sorted(glycan_dict.items(), key=lambda x: monosaccharide_library[x[0]]['mass'])
    
    # Calculate B ions in order
    for i, (sugar, count) in enumerate(sorted_sugars, start=1):
        for j in range(count):
            cumulative += monosaccharide_library[sugar]['mass']
            B_ions[f'B{i}_{j+1}'] = round((cumulative + proton) / charge, 4)

    # --- Calculate Oxonium ions ---
    oxonium_ions = {}
    for sugar in sorted(glycan_dict.keys(), key=lambda x: monosaccharide_library[x]['mass']):
        if glycan_dict[sugar] > 0:
            oxonium_ions[f'ox_{sugar}'] = round(monosaccharide_library[sugar]['mass'] + proton, 4)

    return {
        'b': peptide_ions['b'],
        'y': peptide_ions['y'],
        'c': peptide_ions['c'],
        'z': peptide_ions['z'],
        'immonium': dict(sorted(immonium_ions.items(), key=lambda item: item[1])),
        'internal': dict(sorted(internal_ions.items(), key=lambda item: item[1])),
        'Y': dict(sorted(Y_ions.items(), key=lambda item: item[1])),
        'Y_peptide_loss': dict(sorted(Y_peptide_loss.items(), key=lambda item: item[1])),
        'B': dict(sorted(B_ions.items(), key=lambda item: item[1])),
        'oxonium': dict(sorted(oxonium_ions.items(), key=lambda item: item[1]))
    }

# Add this function to write the results to a CSV file
def write_csv(output_file, data):
    """Writes results to a CSV file."""
    with open(output_file, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["ProteinID", "Site", "Peptide", "Start", "End", "Length", "Sequon", "PredictedMass", "Hydrophobicity", "pI"])
        writer.writeheader()
        writer.writerows(data)

## Experimental -- testing
#from itertools import permutations
def generate_all_y_ions(peptide, glycan_composition, charge=1):
    """
    Generate all possible Y ion series considering different fragmentation paths.

    Parameters:
      peptide (str): Peptide sequence (e.g., "NTSK").
      glycan_composition (str): Glycan composition string (e.g., "HexNAc(5)Hex(5)dHex(1)NeuAc(2)").
      charge (int): The charge state (default is 1).

    Returns:
      dict: A dictionary where keys are fragmentation paths and values are Y ion series.
    """
    # Constants
    proton = 1.007276
    water  = 18.010565

    # Compute the neutral mass of the peptide
    peptide_mass = sum(amino_acid_masses[aa] for aa in peptide) + water

    # Parse glycan composition
    glycan_dict = {}
    for part in glycan_composition.split(')'):
        if part:
            sugar, count = part.split('(')
            glycan_dict[sugar] = int(count)

    # Convert glycan composition into a list of sugars
    glycan_list = []
    for sugar, count in glycan_dict.items():
        glycan_list.extend([sugar] * count)  # Expand each sugar into individual occurrences

    # Generate all possible fragmentation paths
    unique_permutations = set(permutations(glycan_list))  # Unique orders only

    # Compute Y ions for each fragmentation path
    all_Y_ions = {}

    for perm in unique_permutations:
        path_name = " -> ".join(perm)  # Name this path

        current_mass = peptide_mass  # Start with peptide only
        y_series = {'Y0': round((current_mass + proton) / charge, 4)}

        for i, sugar in enumerate(perm, start=1):
            current_mass += monosaccharide_library[sugar]['mass']
            y_series[f'Y{i}'] = round((current_mass + proton) / charge, 4)

        all_Y_ions[path_name] = y_series

    return all_Y_ions

def scan_glycosites(file, glycosylation_type):
    """Scans protein sequences for glycosylation sites and returns site information.
    
    Parameters:
        file (str): Path to the FASTA file
        glycosylation_type (str): Type of glycosylation to scan for (N, O, or C)
        
    Returns:
        list: List of dictionaries containing glycosite information
    """
    glycosites = []
    
    # Regular expression to capture OS and OX from the header
    os_ox_pattern = r"OS=([^\s]+(?: [^\s]+)*)\s+OX=(\d+)\s+GN=([^\s]+)\s+PE=(\d+)\s+SV=(\d+)"
    
    # Process each record in the FASTA file
    for record in SeqIO.parse(file, "fasta"):
        protein_id = record.id
        sequence = str(record.seq)
        header = record.description
        
        # Get metadata from header
        match = re.search(os_ox_pattern, header)
        if match:
            os_value = match.group(1)
            ox_value = match.group(2)
            gn_value = match.group(3)
        else:
            os_value = ""
            ox_value = ""
            gn_value = ""
        
        # Find glycosylation sites
        glyco_sequon = re.compile(glycosylation[glycosylation_type])
        for match in glyco_sequon.finditer(sequence):
            site = match.start() + 1  # 1-based indexing
            sequon = sequence[site-1:site+2]  # Extract sequon
            
            # Get the 21-amino acid sequence centered on the glycosite
            # Calculate start and end positions (0-based)
            start_pos = site - 11  # 10 residues before the site
            end_pos = site + 9     # 10 residues after the site
            
            # Create padded sequence
            seq21 = ""
            if start_pos < 0:
                # Add padding at the start
                seq21 += "Z" * abs(start_pos)
                seq21 += sequence[0:end_pos+1]
            elif end_pos >= len(sequence):
                # Add padding at the end
                seq21 += sequence[start_pos:]
                seq21 += "Z" * (end_pos - len(sequence) + 1)
            else:
                # No padding needed
                seq21 = sequence[start_pos:end_pos+1]
            
            # Ensure the sequence is exactly 21 amino acids
            if len(seq21) != 21:
                if len(seq21) < 21:
                    seq21 += "Z" * (21 - len(seq21))
                else:
                    seq21 = seq21[:21]
            
            glycosites.append({
                "ProteinID": protein_id,
                "Site": site,
                "Sequon": sequon,
                "Species": os_value,
                "TaxonID": ox_value,
                "GeneName": gn_value,
                "seq21": seq21
            })
    
    return glycosites

def write_glycosites_csv(output_file, data):
    """Writes glycosite results to a CSV file."""
    with open(output_file, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["ProteinID", "Site", "Sequon", "Species", "TaxonID", "GeneName", "seq21"])
        writer.writeheader()
        writer.writerows(data)

# Add this to the main function to set up logging
def main():
    
    # PARAMETER SETUP ARGUMENT PARSER
    
    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Glycopeptide Finder")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file. Can be found in the test_proteomes folder.")
    parser.add_argument("-g", "--glycosylation", default="N", help="Glycosylation type (N, O, or C). Default is N. Large file sizes may result from selecting O or C.")
    parser.add_argument("-o", "--output", help="Output CSV file prefix. Default output directory for files is 'data/digested_glycopeptide_library'.")
    parser.add_argument("-p", "--protease", default="trypsin", help="Protease to use for cleavage ('all' for all proteases). Default is trypsin. Proteases: trypsin, chymotrypsin, glu-c, lys-c, arg-c, pepsin, asp-n, proteinase-k.")
    parser.add_argument("-c", "--missed_cleavages", type=int, default=0, help="Number of missed cleavages allowed. Default is 0.")
    parser.add_argument("-m", "--peptide_max_length", type=int, default=25, help="Max peptide length from digestion(default is 25).")
    parser.add_argument("--min-length", type=int, default=5, help="Minimum peptide length (default: 5).")
    parser.add_argument("-y", "--glycan", nargs='?', const='default', help="Path to glycan file (CSV) or use default glycan library. If -y is used without a value, default library will be used.")
    parser.add_argument("-l", "--log", help="Provide log file name. (suggestion: -l log.txt)")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose output.")
    parser.add_argument("-z", "--charge", type=int, default=3, help="Maximum charge state (default: 3).")
    parser.add_argument("--ion-series", action="store_true", help="Compute ion series for peptides and glycopeptides.")

    # Parse arguments
    args = parser.parse_args()

    # Set up logging if log file is provided
    if args.log:
        setup_logging(args.log)

    # Process the input file
    input_file = args.input
    # Extract just the filename without path and extension
    base_filename = os.path.splitext(os.path.basename(input_file))[0]
    missed_cleavages = args.missed_cleavages
    peptide_max_length = args.peptide_max_length
    glycosylation_type = args.glycosylation
    charge_state = args.charge
    glycan_library = args.glycan

    # Get the project root directory (one level up from src)
    project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    # Create output directories
    output_dir = os.path.join(project_root, "data", "digested_glycopeptide_library")
    peptide_output_dir = os.path.join(project_root, "data", "digested_peptide_library")
    glycosite_output_dir = os.path.join(project_root, "data", "glycosite_library")
    
    for directory in [output_dir, peptide_output_dir, glycosite_output_dir]:
        os.makedirs(directory, exist_ok=True)

    # Scan for glycosites and write to file
    glycosites = scan_glycosites(input_file, glycosylation_type)
    glycosite_output_file = f"{glycosite_output_dir}/{base_filename}_{glycosylation_type}-glycosites.csv"
    write_glycosites_csv(glycosite_output_file, glycosites)

    # Create glycosite summary
    glycosite_summary = {}
    for site in glycosites:
        protein_id = site["ProteinID"]
        if protein_id not in glycosite_summary:
            glycosite_summary[protein_id] = 1
        else:
            glycosite_summary[protein_id] += 1
    
    # Write glycosite summary
    summary_file = f"{glycosite_output_dir}/{base_filename}_{glycosylation_type}-glycosites_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Glycosite Summary for {base_filename}\n")
        f.write(f"Glycosylation Type: {glycosylation_type}\n")
        f.write(f"Total Glycosites Found: {len(glycosites)}\n")
        f.write(f"Total Proteins with Glycosites: {len(glycosite_summary)}\n\n")
        f.write("Glycosites per Protein:\n")
        for protein, count in glycosite_summary.items():
            f.write(f"{protein}: {count} glycosites\n")

    if args.verbose:
        print(f"Found {len(glycosites)} glycosites. Results written to {glycosite_output_file}")
        print(f"Glycosite summary written to {summary_file}")
    if args.log:
        logging.info(f"Found {len(glycosites)} glycosites. Results written to {glycosite_output_file}")
        logging.info(f"Glycosite summary written to {summary_file}")

    # Process the input file with the selected protease(s)
    all_results = []

    # If "all" is selected, process all proteases
    if args.protease.lower() == "all":
        selected_proteases = list(proteases.keys())  # All available proteases
    elif args.protease.lower() in proteases:
        selected_proteases = [args.protease.lower()]  # Single protease
    else:
        if args.log:  # Log the error
            logging.error(f"Protease {args.protease} is not supported. Supported proteases: {', '.join(proteases.keys())}")
        return

    # WORKFLOW STARTS HERE

    # Process the input file with the selected protease(s)
    for protease in selected_proteases:
        
        # Log the start of the process
        print(f"Processing {input_file} with protease {protease} and {missed_cleavages} missed cleavages...")
        if args.log:
            logging.info(f"Processing {input_file} with protease {protease} and {missed_cleavages} missed cleavages...")
        
        # Process the input file and return pandas DataFrame
        results_df = process_fasta(input_file, protease, missed_cleavages, glycosylation_type)

        # Flatten results_df to extract peptides
        digest_peptide_library = pd.DataFrame([pep for sublist in results_df["Peptides"].tolist() for pep in sublist])
        
        # Set name of first column to "Peptide"
        digest_peptide_library.columns = ["Peptide"]

        # Add other columns to the DataFrame
        digest_peptide_library["ProteinID"] = results_df["ProteinID"].tolist()[0]

        # Remove empty peptide entries in digest_peptide_library
        digest_peptide_library = digest_peptide_library[digest_peptide_library["Peptide"].str.strip().astype(bool)]

        # Compute peptide mass, hydrophobicity, pI, and m/z values per row
        digest_peptide_library["PredictedMass"] = digest_peptide_library["Peptide"].apply(calculate_peptide_mass)
        digest_peptide_library["Hydrophobicity"] = digest_peptide_library["Peptide"].apply(predict_hydrophobicity)
        digest_peptide_library["pI"] = digest_peptide_library["Peptide"].apply(calculate_pI)
        digest_peptide_library["PredictedMass"] = pd.to_numeric(digest_peptide_library["PredictedMass"], errors='coerce')
        digest_peptide_library = digest_peptide_library.dropna(subset=["PredictedMass"])

        # Write results_df to a CSV file in peptide_library folder
        peptide_output_file = f"{peptide_output_dir}/{base_filename}_{protease}_digested_mc{missed_cleavages}_peptides.csv"
        digest_peptide_library.to_csv(peptide_output_file, index=False)

        if args.verbose:
            print(f"Peptide digestion results written to {peptide_output_file}")
        if args.log:
            logging.info(f"Peptide digestion results written to {peptide_output_file}")

        # Extract the protein sequence and find glycopeptides
        all_glycopeptides = []
        for _, row in results_df.iterrows():
            protein_id = row["ProteinID"]
            sequence = row["Sequence"]

            # Find glycopeptides
            x_glycopeptides = find_glycopeptides(pd.DataFrame([row]), glycosylation_type)
            
            # Process each glycopeptide
            for peptide, site in x_glycopeptides:
                mass = calculate_peptide_mass(peptide)
                hydrophobicity = predict_hydrophobicity(peptide)
                pI = calculate_pI(peptide)
                start_pos = sequence.find(peptide) + 1  # 1-based indexing
                end_pos = start_pos + len(peptide) - 1
                all_glycopeptides.append({
                    "ProteinID": protein_id,
                    "Site": int(site),
                    "Peptide": peptide,
                    "Start": int(start_pos),
                    "End": int(end_pos),
                    "Length": len(peptide),
                    "Sequon": sequence[site - 1:site + 2],  # Extract the sequon amino acid sequence + 1 flanking residue
                    "PredictedMass": mass,
                    "Hydrophobicity": hydrophobicity,
                    "pI": pI,
                })
        
        # Convert to pandas dataframe
        all_glycopeptides = pd.DataFrame(all_glycopeptides)
        
        # Filter by peptide length
        if not all_glycopeptides.empty and 'Peptide' in all_glycopeptides.columns:
            all_glycopeptides = all_glycopeptides.loc[
                (all_glycopeptides['Peptide'].str.len() <= peptide_max_length) &
                (all_glycopeptides['Peptide'].str.len() >= args.min_length)
            ]
        else:
            if args.verbose:
                print("No glycopeptides found; skipping peptide length filtering.")

        # Write the glycopeptide results
        output_file = args.output or f"{output_dir}/{base_filename}_{protease}_digested_mc{missed_cleavages}_z{charge_state}_{glycosylation_type}-glycopeptides.csv"
        all_glycopeptides.to_csv(output_file, index=False)

        if args.verbose:
            print(f"Found {len(all_glycopeptides)} glycopeptides. Results written to {output_file}")
        if args.log:
            logging.info(f"Found {len(all_glycopeptides)} glycopeptides. Results written to {output_file}")

        # Process glycopeptides with glycans if -y flag is used
        if args.glycan is not None:
            try:
                # Use provided glycan file or default based on glycosylation type
                if args.glycan == "default":
                    if glycosylation_type == "N":
                        glycans = default_n_glycan_library
                    elif glycosylation_type == "O":
                        glycans = default_o_glycan_library
                    elif glycosylation_type == "C":
                        glycans = default_c_glycan_library
                    else:
                        raise ValueError(f"No default glycan library for glycosylation type: {glycosylation_type}")
                else:
                    glycans = pd.read_csv(args.glycan)
                
                # Process glycopeptides and compute m/z values
                glycopeptide_results = process_glycopeptides(output_file, glycans, charge_state)

                # Add charge_state from input columns to the DataFrame
                glycopeptide_results["Charge"] = charge_state

                # Compute IonSeries if requested
                if args.ion_series and not glycopeptide_results.empty:
                    glycopeptide_results["IonSeries"] = glycopeptide_results.apply(
                        lambda row: calculate_n_glycopeptide_ions(row["Peptide"], row["Composition"], charge=1),
                        axis=1
                    )
                    output_suffix = "_with_glycans_and_ions"
                else:
                    output_suffix = "_with_glycans"

                # Write the results to a new CSV file
                glycopeptide_output_file = f"{output_dir}/{base_filename}_{protease}_digested_mc{missed_cleavages}_z{charge_state}_{glycosylation_type}-glycopeptides.csv"
                glycopeptide_results.to_csv(glycopeptide_output_file, index=False)

                if args.verbose:
                    print(f"Glycopeptide results written to {glycopeptide_output_file}")
                if args.log:
                    logging.info(f"Glycopeptide results written to {glycopeptide_output_file}")

            except Exception as e:
                if args.log:
                    logging.error(f"Error processing glycan library: {str(e)}")
                if args.verbose:
                    print(f"Error processing glycan library: {str(e)}")
        # If only ion series is requested (no glycans)
        elif args.ion_series and not all_glycopeptides.empty:
            all_glycopeptides["IonSeries"] = all_glycopeptides.apply(
                lambda row: calculate_peptide_ions(row["Peptide"], charge=1),
                axis=1
            )
            # Write results with ion series
            ion_series_output = f"{output_dir}/{base_filename}_{protease}_digested_mc{missed_cleavages}_z{charge_state}_{glycosylation_type}-glycopeptides_with_ions.csv"
            all_glycopeptides.to_csv(ion_series_output, index=False)

# main function
if __name__ == "__main__":
    main()
