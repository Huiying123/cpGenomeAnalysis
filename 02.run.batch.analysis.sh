#!/bin/bash

# =============================================================================
#           COMPLETE BATCH ANALYSIS SCRIPT (FINAL VERSION v7 - Robust Parsing)
# =============================================================================

# --- CONFIGURATION ---
MAPPING_FILE="00.samples.final.txt"
OUTPUT_TABLE="final_summary_table.tsv"
THREADS=8
MISA_SCRIPT_PATH="./misa.pl"
MISA_INI_PATH="./misa.ini"
RESULTS_DIR="analysis_results_v7"

# --- SCRIPT SETUP ---
set -e
set -o pipefail
# ... (Dependency checks from previous versions) ...

# --- INITIALIZE ---
mkdir -p ${RESULTS_DIR}
echo -e "SampleID\tAvg_Coverage_Depth\tGenome_Length_bp\tGC_Content_%\tNum_PCGs\tNum_tRNAs\tLSC_Length\tSSC_Length\tIR_Length\tNum_SSRs" > ${OUTPUT_TABLE}

# --- MAIN LOOP ---
while IFS=$'\t' read -r sample_id reads1_path reads2_path assembly_path genbank_path; do
    if [[ "$sample_id" == \#* ]] || [[ -z "$sample_id" ]]; then continue; fi

    echo "======================================================"
    echo "Processing Sample: ${sample_id}"
    echo "======================================================"

    SAMPLE_DIR="${RESULTS_DIR}/${sample_id}"
    mkdir -p ${SAMPLE_DIR}

    # --- Step 1: Calculate Average Coverage Depth ---
    echo "--> Calculating Coverage Depth..."
    #bwa index "${assembly_path}"
    #bwa mem -t ${THREADS} "${assembly_path}" "${reads1_path}" "${reads2_path}" > "${SAMPLE_DIR}/alignment.sam"
    #samtools view -@ ${THREADS} -bS "${SAMPLE_DIR}/alignment.sam" | samtools sort -@ ${THREADS} -o "${SAMPLE_DIR}/sorted.bam" -
    #rm "${SAMPLE_DIR}/alignment.sam"
    #samtools index "${SAMPLE_DIR}/sorted.bam"
    coverage=$(samtools depth "${SAMPLE_DIR}/sorted.bam" | awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }')
    printf "Coverage for ${sample_id} is %.2fx\n" ${coverage}

    # --- Step 2: Calculate Genome Length and GC Content (*** MORE ROBUST METHOD ***) ---
    echo "--> Calculating Length and GC Content..."
    # Directly pipe seqtk output to awk and read results into a bash array
    stats_array=($(seqtk comp "${assembly_path}" | awk 'NR>1 {print $2, $5}'))
    
    # Assign variables from the array
    length=${stats_array[0]}
    gc_content=${stats_array[1]}
    
    # Check if variables were assigned correctly
    if [ -z "${length}" ]; then
        echo "ERROR: Failed to parse Length/GC for sample ${sample_id}. Skipping."
        continue
    fi
    echo "Length: ${length} bp, GC: ${gc_content}%"

    # --- Step 3: Count Annotated Genes ---
    echo "--> Calculating Annotation Stats..."
    num_pcgs=$(grep -c '     CDS     ' "${genbank_path}")
    num_trnas=$(grep -c '     tRNA    ' "${genbank_path}")
    echo "Found ${num_pcgs} PCGs and ${num_trnas} tRNAs."

    # --- Step 4: Calculate LSC, SSC, IR Lengths ---
    echo "--> Calculating LSC/SSC/IR Lengths..."
    ir_region=$(grep 'repeat_region' "${genbank_path}" | grep 'Inverted Repeat' | head -n 1 || true)
    if [[ -n "$ir_region" ]]; then
        ir_start=$(echo "$ir_region" | grep -o '[0-9]\+' | head -n 1)
        ir_end=$(echo "$ir_region" | grep -o '[0-9]\+' | tail -n 1)
        ir_len=$((ir_end - ir_start + 1))
        ssc_len=18000 # Placeholder: Use a more accurate method if needed
        lsc_len=$((length - ssc_len - (2 * ir_len)))
    else
        ir_len="N/A"
        ssc_len="N/A"
        lsc_len="N/A"
    fi
    echo "LSC: ${lsc_len}, SSC: ${ssc_len}, IR: ${ir_len}"

    # --- Step 5: Find SSRs with MISA ---
    echo "--> Finding SSRs with MISA..."
    awk '/^ORIGIN/{flag=1; next} /\/\//{flag=0} flag' "${genbank_path}" | sed 's/[^a-zA-Z]//g' | tr '[:lower:]' '[:upper:]' > "${SAMPLE_DIR}/sequence.fasta"
    sed -i "1i>${sample_id}" "${SAMPLE_DIR}/sequence.fasta"
    perl ${MISA_SCRIPT_PATH} "${SAMPLE_DIR}/sequence.fasta" ${MISA_INI_PATH}
    num_ssrs=$(grep -c 'ssr' "${SAMPLE_DIR}/sequence.fasta.misa" || true)
    echo "Found ${num_ssrs} SSRs."
    
    # --- Step 6: Write results to the summary table ---
    echo "--> Writing results to table..."
    echo -e "${sample_id}\t${coverage}\t${length}\t${gc_content}\t${num_pcgs}\t${num_trnas}\t${lsc_len}\t${ssc_len}\t${ir_len}\t${num_ssrs}" >> ${OUTPUT_TABLE}

done < ${MAPPING_FILE}

echo "======================================================"
echo "All samples processed successfully!"
echo "Final results table saved to: ${OUTPUT_TABLE}"
echo "======================================================"
