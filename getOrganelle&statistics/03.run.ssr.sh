#!/bin/bash
# FINAL SCRIPT (v13 - Replaced mreps with a self-contained Python script)

# --- CONFIGURATION ---
MAPPING_FILE="00.samples.final3.txt"
OUTPUT_TABLE="ssr_summary_final.tsv"
RESULTS_DIR="ssr_analysis_results"

# --- SCRIPT SETUP ---
set -e
set -o pipefail

# --- DEPENDENCY CHECK ---
echo "--> Checking for required software..."
if ! command -v "python" &> /dev/null; then
    echo "Error: Python is not found."
    exit 1
fi
echo "All required software found."

# --- INITIALIZE ---
mkdir -p ${RESULTS_DIR}
echo -e "SampleID\tNum_SSRs" > ${OUTPUT_TABLE}

# --- MAIN LOOP ---
while read -r sample_id assembly_path; do
    sample_id=$(echo "${sample_id}" | xargs)
    assembly_path=$(echo "${assembly_path}" | xargs)
    if [[ "$sample_id" == \#* ]] || [[ -z "$sample_id" ]]; then continue; fi

    echo "======================================================"
    echo "Processing Sample: ${sample_id}"
    echo "======================================================"

    if [ ! -f "${assembly_path}" ]; then
        echo "--> ERROR: FASTA file not found: '${assembly_path}'. Skipping."
        continue
    fi

    # --- Find SSRs with our new Python script ---
    echo "--> Finding SSRs with Python script..."
    
    # Call the python script and pass the assembly path as an argument
    num_ssrs=$(python find_ssrs.py "${assembly_path}")
    
    echo "Found ${num_ssrs} SSRs."
    
    echo "--> Writing results to table..."
    echo -e "${sample_id}\t${num_ssrs}" >> ${OUTPUT_TABLE}

done < <(tr -d '\r' < "${MAPPING_FILE}")

echo "======================================================"
echo "All samples processed successfully!"
echo "Final results table saved to: ${OUTPUT_TABLE}"
echo "======================================================"
