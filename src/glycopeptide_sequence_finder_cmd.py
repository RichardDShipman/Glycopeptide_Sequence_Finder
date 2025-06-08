"""
glycopeptide_sequence_finder_refactored.py

This script processes a FASTA file, or all FASTA files in a directory, to 
identify glycopeptides based on protease cleavage rules and glycosylation 
sequons. It predicts their masses, hydrophobicity, and pI, and can calculate 
theoretical ion series. The results are written to CSV files and a final
summary report is generated.

This refactored version encapsulates logic within classes for better organization
and maintainability, while remaining a single script.

Usage:
    python glycopeptide_sequence_finder_refactored.py -i <input_fasta_file_or_directory> -p <protease>

Author:
    Richard Shipman -- 2025
Refactored with:
    Gemini 2.5 Pro, ChatGPT 4o
"""
import argparse
import csv
import logging
import re
from datetime import datetime
from itertools import chain, combinations
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO


class GlycoPipeConfig:
    """
    Configuration class to hold all constants and biological data.
    """

    PROTEASES: Dict[str, Tuple[str, Optional[str]]] = {
        "trypsin": ("[KR]", "P"),
        "chymotrypsin": ("[FLWY]", "P"),
        "glu-c": ("[DE]", "P"),
        "lys-c": ("K", "P"),
        "arg-c": ("R", "P"),
        "pepsin": ("[FWY]", None),
        "asp-n": ("D", None),
        "proteinase-k": ("[AFILVWY]", None),
        "operator": ("[ST]", None),
        "thermolysin": ("[ALIVFM]", "[DE]"),
    }

    GLYCOSYLATION: Dict[str, str] = {
        "N": "N[^P][ST]",
        "O": "[ST]",
        "C": "W..[WCF]",
        "C-Mannose-1": "W..W..W",
        "C-Mannose-2": "W..W..W..C",
        "C-Mannose-3": "W..[WC]",
        "N-Glycan": "N[^P][STC]",
        "O-Fucose-1": "C....[ST]C",
        "O-Fucose-2": "C..[ST]C",
        "O-Fucose-3": "C..[ST]C..G",
        "O-GlcNAc": "C....[ST]G..C",
        "O-Glucose-1": "C.NT.GS[FY].C",
        "O-Glucose-2": "C.S.[PA]C",
        "type-mutation-c": "[^N][^P][ST]",
        "type-mutation-b": "NP[ST]",
        "type-mutation-a": "N[^P][^S^T]",
        "type-extended-mutation-a": "[^N][^P][STC]",
        "type-extended-mutation-b": "NP[STC]",
        "type-extended-mutation-c": "N[^P][^S^T^C]",
    }

    AMINO_ACID_MASSES: Dict[str, float] = {
        'A': 71.03711, 'R': 156.10111, 'N': 114.04293, 'D': 115.02694,
        'C': 103.00919, 'Q': 128.05858, 'E': 129.04259, 'G': 57.02146,
        'H': 137.05891, 'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
        'M': 131.04049, 'F': 147.06841, 'P': 97.05276, 'S': 87.03203,
        'T': 101.04768, 'W': 186.07931, 'Y': 163.06333, 'V': 99.06841
    }
    
    WATER_MASS = 18.010565
    PROTON_MASS = 1.007276
    NH3_MASS = 17.0265

    HYDROPHOBICITY_VALUES: Dict[str, float] = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5,
        'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
        'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
        'Y': -1.3, 'V': 4.2
    }

    PKA_VALUES: Dict[str, float] = {
        'C': 8.18, 'D': 3.65, 'E': 4.25, 'H': 6.00, 'K': 10.53, 'R': 12.48,
        'Y': 10.07, 'N': 3.22, 'Q': 3.22, 'S': 3.70, 'T': 3.70, 'W': 10.07
    }
    
    TERMINAL_PKA = {'N_term': 9.6, 'C_term': 2.3}

    DEFAULT_N_GLYCAN_LIBRARY = pd.DataFrame([
        {"glytoucan_ac": "G22768VO", "byonic": "HexNAc(2)Hex(3) % 1216.422863", "composition": "HexNAc(2)Hex(3)", "mass": 1216.422863, "shorthand_glycan": "N2H3"},
    ])

    DEFAULT_O_GLYCAN_LIBRARY = pd.DataFrame([
        {"glytoucan_ac": "G14843DJ", "byonic": "HexNAc(1) % 221.089937305", "composition": "HexNAc(1)", "mass": 221.089937305, "shorthand_glycan": "N1"},
    ])

    DEFAULT_C_GLYCAN_LIBRARY = pd.DataFrame([
        {"glytoucan_ac": "G81399MY", "byonic": "Hex(1) % 180.0633882", "composition": "Hex(1)", "mass": 180.0633882, "shorthand_glycan": "H1"},
    ])

    MONOSACCHARIDE_LIBRARY: Dict[str, Dict[str, Any]] = {
        "Hex": {"mass": 162.0528, "formula": "C6H10O5", "symbol": "H"},
        "HexNAc": {"mass": 203.0794, "formula": "C8H13NO5", "symbol": "N"},
        "dHex": {"mass": 146.0579, "formula": "C6H12O5", "symbol": "dH"},
        "NeuAc": {"mass": 291.0954, "formula": "C11H17NO8", "symbol": "S"},
    }
    
    IMMONIUM_MASSES: Dict[str, float] = {
        'A': 44.05003, 'R': 129.11400, 'N': 87.05584, 'D': 88.03986, 'C': 76.02210,
        'E': 102.05550, 'Q': 101.07150, 'G': 30.03438, 'H': 110.07180, 'I': 86.09698,
        'L': 86.09698, 'K': 101.10790, 'M': 104.05340, 'F': 120.08130, 'P': 70.06568,
        'S': 60.04494, 'T': 74.06059, 'W': 159.09220, 'Y': 136.07620, 'V': 72.08133
    }

class GlycopeptideProcessor:
    """
    Orchestrates the glycopeptide finding pipeline from a FASTA file.
    """
    def __init__(self, args: argparse.Namespace):
        self.args = args
        self.config = GlycoPipeConfig()
        self.project_root = Path(__file__).resolve().parent.parent
        self.output_dir: Path
        self.peptide_output_dir: Path
        self.glycosite_output_dir: Path
        self.stats: Dict[str, Any] = {'file_name': Path(self.args.input).name, 'protease_stats': {}}

    def run(self) -> Dict[str, Any]:
        """Executes the full pipeline for a single file and returns stats."""
        self._setup_directories()
        
        base_filename = Path(self.args.input).stem
        self.stats['protein_count'] = self._count_proteins_in_fasta()

        glycosites = self._scan_and_write_glycosites(base_filename)
        self.stats['total_glycosites'] = len(glycosites)
        
        selected_proteases = self._get_selected_proteases()
        if not selected_proteases:
            logging.error(f"Protease {self.args.protease} is not supported.")
            return self.stats

        for protease in selected_proteases:
            protease_stats: Dict[str, Any] = {}
            logging.info(f"Processing with protease '{protease}'...")
            print(f"Processing '{base_filename}' with protease '{protease}'...")
            
            digested_proteins_df = self._process_fasta_to_peptides(protease)
            
            glycopeptide_candidates_df = self._find_glycopeptides_from_digest(digested_proteins_df)
            protease_stats['initial_glycopeptides'] = len(glycopeptide_candidates_df)
            
            glycopeptide_candidates_df = self._filter_peptides_by_length(glycopeptide_candidates_df)
            protease_stats['filtered_glycopeptides'] = len(glycopeptide_candidates_df)
            
            if glycopeptide_candidates_df.empty:
                logging.warning(f"No glycopeptides found for protease {protease} after filtering. Skipping.")
                self.stats['protease_stats'][protease] = protease_stats
                continue

            peptide_output_file = self.peptide_output_dir / f"{base_filename}_{protease}_mc{self.args.missed_cleavages}_peptides.csv"
            glycopeptide_candidates_df.to_csv(peptide_output_file, index=False)
            logging.info(f"Peptide library saved to {peptide_output_file}")

            final_df = self._process_glycans_and_ions(glycopeptide_candidates_df)
            protease_stats['final_glycopeptide_variants'] = len(final_df)

            suffix = "glycopeptides_with_glycans_and_ions" if self.args.ion_series and self.args.glycan else \
                     "glycopeptides_with_glycans" if self.args.glycan else \
                     "glycopeptides_with_ions" if self.args.ion_series else "glycopeptides"

            output_file = self.output_dir / f"{base_filename}_{protease}_mc{self.args.missed_cleavages}_z{self.args.charge}_{self.args.glycosylation}_{suffix}.csv"
            final_df.to_csv(output_file, index=False)
            logging.info(f"Final results for protease '{protease}' saved to {output_file}")
            print(f"Final results for '{protease}' saved to {output_file}")
            
            self.stats['protease_stats'][protease] = protease_stats
        
        return self.stats

    def _count_proteins_in_fasta(self) -> int:
        """Counts the number of records in the input FASTA file."""
        try:
            with open(self.args.input, "r") as f:
                return sum(1 for line in f if line.startswith('>'))
        except Exception as e:
            logging.error(f"Could not count proteins in {self.args.input}: {e}")
            return 0

    def _setup_directories(self):
        """Creates required output directories."""
        self.output_dir = self.project_root / "data" / "digested_glycopeptide_library"
        self.peptide_output_dir = self.project_root / "data" / "digested_peptide_library"
        self.glycosite_output_dir = self.project_root / "data" / "glycosite_library"
        self.reports_dir = self.project_root / "reports"
        
        for directory in [self.output_dir, self.peptide_output_dir, self.glycosite_output_dir, self.reports_dir]:
            directory.mkdir(parents=True, exist_ok=True)
        logging.info("Output directories created/verified.")
        
    def _get_selected_proteases(self) -> List[str]:
        """Returns a list of proteases to use based on user arguments."""
        protease_arg = self.args.protease.lower()
        if protease_arg == "all":
            return list(self.config.PROTEASES.keys())
        if protease_arg in self.config.PROTEASES:
            return [protease_arg]
        return []

    def _cleave_sequence(self, sequence: str, protease_name: str, missed_cleavages: int) -> List[str]:
        """Cleaves a sequence based on protease rules."""
        cleavage_pattern, exclusion = self.config.PROTEASES[protease_name]

        regex = rf"(?<={cleavage_pattern})(?!{exclusion})" if protease_name != "asp-n" else rf"(?={cleavage_pattern})"
        
        fragments = re.split(regex, sequence)
        unique_peptides = set()

        for i in range(len(fragments)):
            for j in range(i + 1, min(i + 2 + missed_cleavages, len(fragments) + 1)):
                peptide = "".join(fragments[i:j])
                if len(peptide) >= self.args.min_length:
                    unique_peptides.add(peptide)

        return sorted(list(unique_peptides))

    def _process_fasta_to_peptides(self, protease: str) -> pd.DataFrame:
        """Processes a FASTA file to generate a DataFrame of digested proteins."""
        results = []
        os_ox_pattern = re.compile(r"OS=([^\s]+(?: [^\s]+)*)\s+OX=(\d+)\s+GN=([^\s]+)\s+PE=(\d+)\s+SV=(\d+)")

        for record in SeqIO.parse(self.args.input, "fasta"):
            header = record.description
            match = os_ox_pattern.search(header)
            
            peptides = self._cleave_sequence(str(record.seq), protease, self.args.missed_cleavages)
            
            results.append({
                "ProteinID": record.id,
                "Peptides": peptides,
                "Sequence": str(record.seq),
                "Species": match.group(1) if match else "",
                "TaxonID": match.group(2) if match else "",
                "GeneName": match.group(3) if match else "",
            })
        return pd.DataFrame(results)

    def _find_glycopeptides_from_digest(self, digested_df: pd.DataFrame) -> pd.DataFrame:
        """Identifies peptides with glycosylation sequons from a digest."""
        glyco_sequon = re.compile(self.config.GLYCOSYLATION[self.args.glycosylation])
        all_glycopeptides = []

        for _, row in digested_df.iterrows():
            for peptide in row['Peptides']:
                for match in glyco_sequon.finditer(peptide):
                    site_in_protein = row['Sequence'].find(peptide) + match.start() + 1
                    all_glycopeptides.append({
                        "ProteinID": row['ProteinID'],
                        "Site": site_in_protein,
                        "Peptide": peptide,
                        "Start": row['Sequence'].find(peptide) + 1,
                        "End": row['Sequence'].find(peptide) + len(peptide),
                        "Length": len(peptide),
                        "Sequon": peptide[match.start():match.end()],
                    })

        if not all_glycopeptides:
            return pd.DataFrame()

        glyco_df = pd.DataFrame(all_glycopeptides).drop_duplicates()
        glyco_df["PredictedMass"] = glyco_df["Peptide"].apply(self._calculate_peptide_mass)
        glyco_df["Hydrophobicity"] = glyco_df["Peptide"].apply(self._predict_hydrophobicity)
        glyco_df["pI"] = glyco_df["Peptide"].apply(self._calculate_pI)
        
        return glyco_df

    def _filter_peptides_by_length(self, df: pd.DataFrame) -> pd.DataFrame:
        """Filters a DataFrame of peptides based on min and max length."""
        if df.empty:
            return df
        return df[df['Peptide'].str.len() <= self.args.peptide_max_length]

    @staticmethod
    def _calculate_peptide_mass(sequence: str) -> Optional[float]:
        """Calculates the mass of a peptide."""
        if any(aa not in GlycoPipeConfig.AMINO_ACID_MASSES for aa in sequence):
            return None
        mass = sum(GlycoPipeConfig.AMINO_ACID_MASSES.get(aa, 0) for aa in sequence) + GlycoPipeConfig.WATER_MASS
        return mass

    @staticmethod
    def _predict_hydrophobicity(peptide_sequence: str) -> float:
        """Predicts hydrophobicity using the Kyte-Doolittle scale."""
        if not peptide_sequence or any(aa not in GlycoPipeConfig.HYDROPHOBICITY_VALUES for aa in peptide_sequence):
            return 0.0
        total_hydrophobicity = sum(GlycoPipeConfig.HYDROPHOBICITY_VALUES.get(aa, 0) for aa in peptide_sequence)
        return round(total_hydrophobicity / len(peptide_sequence), 5)

    @staticmethod
    def _calculate_pI(peptide_sequence: str) -> float:
        """Calculates the isoelectric point (pI) of a peptide."""
        # Check only for amino acids that have pKa values to avoid issues with 'X' etc.
        ionizable_residues = set(GlycoPipeConfig.PKA_VALUES.keys())
        if any(aa not in GlycoPipeConfig.AMINO_ACID_MASSES for aa in peptide_sequence):
             return 0.0

        residue_counts = {aa: peptide_sequence.count(aa) for aa in ionizable_residues}

        def _calculate_net_charge(ph: float) -> float:
            charge = 1 / (1 + 10**(ph - GlycoPipeConfig.TERMINAL_PKA['N_term']))
            charge -= 1 / (1 + 10**(GlycoPipeConfig.TERMINAL_PKA['C_term'] - ph))

            for aa, count in residue_counts.items():
                if count > 0:
                    pka = GlycoPipeConfig.PKA_VALUES[aa]
                    if aa in ['D', 'E', 'Y', 'C']:
                        charge -= count / (1 + 10**(pka - ph))
                    elif aa in ['H', 'K', 'R']:
                        charge += count / (1 + 10**(ph - pka))
            return charge

        low, high = 0.0, 14.0
        while high - low > 0.01:
            mid = (low + high) / 2
            if _calculate_net_charge(mid) > 0:
                low = mid
            else:
                high = mid
        return round((low + high) / 2, 2)

    def _process_glycans_and_ions(self, df: pd.DataFrame) -> pd.DataFrame:
        """Adds glycan information and/or ion series to the DataFrame."""
        if self.args.glycan:
            glycan_lib = self._load_glycan_library()
            df = self._add_glycan_data(df, glycan_lib)
            if self.args.ion_series:
                 df["IonSeries"] = df.apply(
                    lambda row: self._calculate_n_glycopeptide_ions(row["Peptide"], row["composition"], charge=1),
                    axis=1
                )
        elif self.args.ion_series:
            df["IonSeries"] = df["Peptide"].apply(self._calculate_peptide_ions, charge=1)
        return df

    def _load_glycan_library(self) -> pd.DataFrame:
        """Loads the appropriate glycan library."""
        if self.args.glycan == "default":
            g_type = self.args.glycosylation
            if g_type == "N": return self.config.DEFAULT_N_GLYCAN_LIBRARY
            if g_type == "O": return self.config.DEFAULT_O_GLYCAN_LIBRARY
            if g_type == "C": return self.config.DEFAULT_C_GLYCAN_LIBRARY
            raise ValueError(f"No default glycan library for type: {g_type}")
        return pd.read_csv(self.args.glycan)

    def _add_glycan_data(self, peptides_df: pd.DataFrame, glycans_df: pd.DataFrame) -> pd.DataFrame:
        """Combines peptide and glycan data, calculating masses and m/z."""
        peptides_df['key'] = 1
        glycans_df['key'] = 1
        
        glycopeptide_df = pd.merge(peptides_df, glycans_df, on='key').drop('key', axis=1)

        glycopeptide_df['GlycopeptideMass'] = glycopeptide_df['PredictedMass'] + glycopeptide_df['mass']
        glycopeptide_df.rename(columns={'mass': 'GlycanMass'}, inplace=True)

        for z in range(2, self.args.charge + 1):
            glycopeptide_df[f'z{z}'] = (glycopeptide_df['GlycopeptideMass'] + (z * self.config.PROTON_MASS)) / z
        
        return glycopeptide_df

    def _scan_and_write_glycosites(self, base_filename: str) -> List[Dict]:
        """Scans for glycosites and writes summary and detailed files."""
        logging.info(f"Scanning for {self.args.glycosylation}-glycosites in {self.args.input}...")
        glycosites = []
        os_ox_pattern = re.compile(r"OS=([^\s]+(?: [^\s]+)*)\s+OX=(\d+)\s+GN=([^\s]+)")

        for record in SeqIO.parse(self.args.input, "fasta"):
            match = os_ox_pattern.search(record.description)
            glyco_sequon = re.compile(self.config.GLYCOSYLATION[self.args.glycosylation])
            
            for m in glyco_sequon.finditer(str(record.seq)):
                glycosites.append({
                    "ProteinID": record.id,
                    "Site": m.start() + 1,
                    "Sequon": m.group(0),
                    "Species": match.group(1) if match else "",
                    "TaxonID": match.group(2) if match else "",
                    "GeneName": match.group(3) if match else "",
                })
        
        if not glycosites:
            logging.warning("No glycosites found.")
            return []

        glycosites_df = pd.DataFrame(glycosites)
        output_file = self.glycosite_output_dir / f"{base_filename}_{self.args.glycosylation}-glycosites.csv"
        glycosites_df.to_csv(output_file, index=False)
        logging.info(f"Glycosite scan results written to {output_file}")
        return glycosites

    def _calculate_peptide_ions(self, peptide: str, charge: int = 1) -> Dict[str, List[float]]:
        """Calculates theoretical m/z values for b, y, c, z ions."""
        if any(aa not in self.config.AMINO_ACID_MASSES for aa in peptide):
            logging.warning(f"Peptide '{peptide}' contains non-standard amino acids. Skipping ion calculation.")
            return {}
            
        ions: Dict[str, List[float]] = {'b': [], 'y': [], 'c': [], 'z': []}
        
        # b and c ions
        cumulative_mass = 0.0
        for i in range(len(peptide) - 1):
            cumulative_mass += self.config.AMINO_ACID_MASSES[peptide[i]]
            ions['b'].append(round((cumulative_mass + self.config.PROTON_MASS) / charge, 4))
            ions['c'].append(round((cumulative_mass + self.config.NH3_MASS + self.config.PROTON_MASS) / charge, 4))

        # y and z ions
        cumulative_mass = 0.0
        for i in range(len(peptide) - 1, 0, -1):
            cumulative_mass += self.config.AMINO_ACID_MASSES[peptide[i]]
            ions['y'].append(round((cumulative_mass + self.config.WATER_MASS + self.config.PROTON_MASS) / charge, 4))
            ions['z'].append(round((cumulative_mass + self.config.PROTON_MASS - self.config.NH3_MASS) / charge, 4))
        
        for ion_type in ions: ions[ion_type].sort()
        return ions
        
    def _calculate_n_glycopeptide_ions(self, peptide: str, glycan_composition: str, charge: int = 1) -> Dict:
        """
        Calculates theoretical ions for an N-glycopeptide, including peptide fragments (b, y, c, z),
        glycan fragments (B-ions), and peptide+glycan fragments (Y-ions).
        """
        final_ions = self._calculate_peptide_ions(peptide, charge)
        if not final_ions:
            return {}

        try:
            glycan_list = []
            if not isinstance(glycan_composition, str) or not glycan_composition.strip():
                logging.warning(f"Empty or invalid glycan composition for peptide {peptide}. Skipping B/Y ion calculation.")
                final_ions['B'] = []
                final_ions['Y'] = []
                return final_ions
                
            parts = re.findall(r"(\w+)\((\d+)\)", glycan_composition)
            for sugar, count in parts:
                if sugar in self.config.MONOSACCHARIDE_LIBRARY:
                    glycan_list.extend([sugar] * int(count))
                else:
                    logging.warning(f"Unknown monosaccharide '{sugar}' in composition '{glycan_composition}'. Skipping.")

            if not glycan_list:
                final_ions['B'] = []
                final_ions['Y'] = []
                return final_ions

            b_ion_masses = set()
            all_subsets = chain.from_iterable(combinations(glycan_list, r) for r in range(1, len(glycan_list) + 1))
            
            for subset in all_subsets:
                mass = sum(self.config.MONOSACCHARIDE_LIBRARY[s]['mass'] for s in subset)
                b_ion_masses.add(mass)
                
            final_ions['B'] = sorted([round(m / charge, 4) for m in b_ion_masses])

            peptide_mass = self._calculate_peptide_mass(peptide)
            if peptide_mass is None:
                final_ions['Y'] = []
                return final_ions
                
            y_ion_masses = {peptide_mass + self.config.MONOSACCHARIDE_LIBRARY['HexNAc']['mass']}
            y_ion_masses.update(peptide_mass + b_mass for b_mass in b_ion_masses)
            final_ions['Y'] = sorted([round((m + self.config.PROTON_MASS) / charge, 4) for m in y_ion_masses])

        except Exception as e:
            logging.error(f"Could not calculate B/Y ions for peptide '{peptide}' with composition '{glycan_composition}': {e}")
            final_ions['B'], final_ions['Y'] = [], []

        return final_ions

def parse_arguments() -> argparse.Namespace:
    """Sets up and parses command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Glycopeptide Finder: A tool for in-silico digestion and glycopeptide identification.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file or directory containing FASTA files.")
    parser.add_argument("-g", "--glycosylation", default="N", choices=list(GlycoPipeConfig.GLYCOSYLATION.keys()), help="Glycosylation type to search for.")
    parser.add_argument("-p", "--protease", default="trypsin", help=f"Protease for digestion. Use 'all' for all proteases. Choices: {', '.join(GlycoPipeConfig.PROTEASES.keys())}")
    parser.add_argument("-c", "--missed_cleavages", type=int, default=0, help="Number of allowed missed cleavages.")
    parser.add_argument("-m", "--peptide_max_length", type=int, default=25, help="Maximum peptide length from digestion.")
    parser.add_argument("--min-length", type=int, default=5, help="Minimum peptide length.")
    parser.add_argument("-y", "--glycan", nargs='?', const='default', help="Path to a glycan CSV file, or use 'default' for a built-in library.")
    parser.add_argument("-l", "--log", help="Log file name. If not provided, no log file will be created.")
    parser.add_argument("-z", "--charge", type=int, default=3, help="Maximum charge state for m/z calculation.")
    parser.add_argument("--ion-series", action="store_true", help="Compute theoretical ion series for peptides/glycopeptides.")
    
    return parser.parse_args()

def generate_summary_report(all_stats: List[Dict[str, Any]], args: argparse.Namespace):
    """Generates and saves a summary report of the run."""
    report_lines = []
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    report_lines.append("="*80)
    report_lines.append("Glycopeptide Finder - Run Summary Report")
    report_lines.append(f"Run Date: {timestamp}")
    report_lines.append("="*80)

    report_lines.append("\n--- Global Parameters ---")
    report_lines.append(f"Protease(s): {args.protease}")
    report_lines.append(f"Glycosylation Type: {args.glycosylation}")
    report_lines.append(f"Missed Cleavages: {args.missed_cleavages}")
    report_lines.append(f"Peptide Length Range: {args.min_length}-{args.peptide_max_length}")
    report_lines.append(f"Glycan Library: {'Default' if args.glycan == 'default' else args.glycan if args.glycan else 'None'}")
    report_lines.append(f"Ion Series Calculation: {'Enabled' if args.ion_series else 'Disabled'}")

    grouped_stats: Dict[str, List[Any]] = {}
    for stats in all_stats:
        file_name = stats.get('file_name', 'UnknownFile')
        group_key = file_name.split('_')[0] if '_' in file_name else file_name
        if group_key not in grouped_stats:
            grouped_stats[group_key] = []
        grouped_stats[group_key].append(stats)
        
    for group, stats_list in grouped_stats.items():
        report_lines.append(f"\n\n--- Species/Group: {group.capitalize()} ---")
        for stats in stats_list:
            report_lines.append(f"\n  File: {stats['file_name']}")
            report_lines.append(f"  - Total Proteins Processed: {stats.get('protein_count', 'N/A')}")
            report_lines.append(f"  - Total Potential Glycosites Found: {stats.get('total_glycosites', 'N/A')}")
            for protease, p_stats in stats.get('protease_stats', {}).items():
                report_lines.append(f"    - Stats for Protease '{protease}':")
                report_lines.append(f"      - Initial Glycopeptide Candidates: {p_stats.get('initial_glycopeptides', 'N/A')}")
                report_lines.append(f"      - Glycopeptides after Length Filter: {p_stats.get('filtered_glycopeptides', 'N/A')}")
                if args.glycan:
                    report_lines.append(f"      - Final Glycopeptide-Glycan Variants: {p_stats.get('final_glycopeptide_variants', 'N/A')}")

    report_content = "\n".join(report_lines)
    print("\n" + report_content)

    project_root = Path(__file__).resolve().parent.parent
    reports_dir = project_root / "reports"
    reports_dir.mkdir(exist_ok=True)
    report_filename = f"run_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    report_path = reports_dir / report_filename
    
    try:
        report_path.write_text(report_content)
        print(f"\nSummary report saved to: {report_path}")
        logging.info(f"Summary report saved to: {report_path}")
    except Exception as e:
        print(f"Error saving summary report: {e}")
        logging.error(f"Error saving summary report: {e}")


def main():
    """Main entry point for the script."""
    args = parse_arguments()

    if args.log:
        logging.basicConfig(
            filename=args.log, level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s', filemode='w'
        )
    
    input_path = Path(args.input)
    if not input_path.exists():
        msg = f"Error: Input path '{input_path}' does not exist."
        print(msg)
        logging.error(msg)
        return

    fasta_files = []
    if input_path.is_dir():
        print(f"Input is a directory. Searching for FASTA files...")
        extensions = ["*.fasta", "*.fa", "*.fna", "*.faa"]
        for ext in extensions:
            fasta_files.extend(input_path.glob(ext))
        if not fasta_files:
            msg = f"No FASTA files found in directory '{input_path}' with extensions {extensions}"
            print(msg)
            logging.warning(msg)
            return
        print(f"Found {len(fasta_files)} FASTA file(s) to process.")
    elif input_path.is_file():
        fasta_files.append(input_path)
    else:
        msg = f"Error: Input path '{input_path}' is not a valid file or directory."
        print(msg)
        logging.error(msg)
        return
        
    all_run_stats = []
    for fasta_file in fasta_files:
        print(f"\n{'='*20}\nProcessing file: {fasta_file.name}\n{'='*20}")
        logging.info(f"--- Processing file: {fasta_file.name} ---")
        
        file_specific_args = argparse.Namespace(**vars(args))
        file_specific_args.input = str(fasta_file)
        
        pipeline = GlycopeptideProcessor(file_specific_args)
        stats = pipeline.run()
        all_run_stats.append(stats)

    if all_run_stats:
        generate_summary_report(all_run_stats, args)

    print("\nPipeline finished for all files.")
    logging.info("--- Pipeline finished for all files. ---")


if __name__ == "__main__":
    main()
