#!/bin/bash
set -e  # Exit on error
trap 'echo "Error occurred. Cleaning up..."; exit 1' ERR

start_time=$(date +%s)

# Default parameters
input_dir="data/test_proteomes"
protease="trypsin"
glycosylation_type="N"
missed_cleavages=0
charge_state=2
max_peptide_length=25
min_peptide_length=5
cores=4

# Get the project root directory (where this script is located)
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Create output directories if they don't exist
mkdir -p "${PROJECT_ROOT}/data/digested_peptide_library"
mkdir -p "${PROJECT_ROOT}/data/digested_glycopeptide_library"
mkdir -p "${PROJECT_ROOT}/data/logs"

# Generate timestamp for this run
timestamp=$(date +%Y%m%d_%H%M%S)
log_file="${PROJECT_ROOT}/data/logs/batch_run_${timestamp}.log"

# Protease (-p)
#protease="trypsin"  # Cleaves after K or R unless followed by P
#protease="chymotrypsin"  # Cleaves after F, W, L, or Y unless followed by P
#protease="pepsin"  # Cleaves after F, W, or Y
#protease="glu-c"  # Cleaves after E
#protease="lys-c"  # Cleaves after K
#protease="arg-c"  # Cleaves after R
#protease="asp-n"  # Cleaves before Asp (D)
#protease="proteinase-k"  # Cleaves after A, F, I, L, V, W, Y
#protease="operator" # cleaves after glycosylated S or T, for O-glycans
#protease="all" # all proteases

# Glycosylation type (-g)
#glycosylation_type="N"
#glycosylation_type="O"
#glycosylation_type="C"

# welcome message
ascii_glycopeptide1="
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  🔬 F I N D I N G   G L Y C O P E P T I D E   S E Q U E N C E S ! ! ! 🔬
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
                 GLYCO-
      H-N-C-C-O--PEPTIDE--N-C-C-O-H-N-C-C-O-H
          |      SEQUENCE   |         |
          R      FINDER     R         R
         /                  \\          \\
        N-Glycan             O-Glycan   C-Glycan
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  📊 P R E D I C T I N G   G L Y C O P E P T I D E   P R O P E R T I E S ! ! ! 📊
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Welcome message
echo "Welcome to the Batch Run of the Glycopeptide Sequence Finder!"
echo "$ascii_glycopeptide1"
echo "Starting glycopeptide sequence finder..."
echo "Please wait while the glycopeptide sequence finder is running..."

# Define the processing function
process_fasta() {
    local fasta_file="$1"
    echo "Processing $fasta_file..."
    
    # Generate output filename
    output_file="${PROJECT_ROOT}/data/digested_glycopeptide_library/$(basename ${fasta_file} .fasta)_${protease}_digested_mc${missed_cleavages}_z${charge_state}_${glycosylation_type}-glycopeptides.csv"
    
    # Run with error handling
    if python "${PROJECT_ROOT}/src/glycopeptide_sequence_finder_cmd.py" \
        -i "${fasta_file}" \
        -o "${output_file}" \
        -p ${protease} \
        -g ${glycosylation_type} \
        -c ${missed_cleavages} \
        -l "${PROJECT_ROOT}/data/logs/$(basename ${fasta_file} .fasta)_${timestamp}.log" \
        -v \
        -z ${charge_state} \
        -m ${max_peptide_length} \
        --min-length ${min_peptide_length} >> "${log_file}"; then
        echo "Successfully processed ${fasta_file}"
    else
        echo "Error processing ${fasta_file}. Check ${log_file} for details."
    fi
}

# Export the function and variables needed by it
export -f process_fasta
export PROJECT_ROOT protease glycosylation_type missed_cleavages charge_state max_peptide_length min_peptide_length timestamp log_file

# Process each FASTA file in parallel using xargs
find "${PROJECT_ROOT}/${input_dir}" -name "*.fasta" | \
xargs -P ${cores} -I {} bash -c 'process_fasta "{}"'

echo "Digested and tasted the glycoproteome. Yummy! 🍽️"

# calculate elapsed time 
end_time=$(date +%s)
elapsed_time=$((end_time - start_time))
elapsed_time_mins=$(echo "scale=2; $elapsed_time / 60" | bc)
echo "Total elapsed time: $elapsed_time_mins minutes."

# Safely count files with error handling
count_entries() {
    local dir=$1
    if [ -d "$dir" ] && [ "$(ls -A $dir)" ]; then
        tail -n +2 "$dir"/*.csv 2>/dev/null | wc -l
    else
        echo "0"
    fi
}

# Count sequences
peptide_count=$(count_entries "${PROJECT_ROOT}/data/digested_peptide_library")
glycopeptide_count=$(count_entries "${PROJECT_ROOT}/data/digested_glycopeptide_library")

# Create summary report
summary_report="${PROJECT_ROOT}/data/logs/summary_batch_run_${timestamp}.txt"

# Write the summary
{
    echo "==================== Glycopeptide Sequence Finder Summary ===================="
    echo "Batch Run Date: $(date)"
    echo "Total Peptide Sequences Found: $peptide_count"
    echo "Total Glycopeptide Sequences Found: $glycopeptide_count"
    echo "---------------------------------------------------------------------------"
    echo "Parameters Used:"
    echo "  Protease: $protease"
    echo "  Missed Cleavages: $missed_cleavages"
    echo "  Glycosylation Type: $glycosylation_type"
    echo "  Min Peptide Length: $min_peptide_length"
    echo "  Max Peptide Length: $max_peptide_length"
    echo "  Charge State: $charge_state"
    echo "---------------------------------------------------------------------------"
    echo "Processing Time: $elapsed_time_mins minutes"
    echo "Log File: $log_file"
    echo "---------------------------------------------------------------------------"
} > "$summary_report"

echo "Summary saved to: $summary_report"
echo "Log file saved to: $log_file"


