#!/bin/bash

# =============================================================================
#           A BASH SCRIPT FOR BATCH PROCESSING OF 64 GENOME SAMPLES
# =============================================================================
#
# PURPOSE:
#   - Calculates key statistics for a list of genome samples.
#   - Generates a summary table suitable for publication.
#
# REQUIREMENTS:
#   - A mapping file (samples_map.txt) with 5 columns:
#     SampleID, Path_Read1, Path_Read2, Path_Assembly, Path_GenBank
#   - Conda environment with bwa, samtools, seqtk, and busco installed.
#
# USAGE:
#   - Update the CONFIGURATION section below.
#   - Run in terminal: ./run_full_analysis.sh
#
# =============================================================================


# --- CONFIGURATION ---
# (请根据您的系统和数据进行修改)

# 输入的样本映射文件
MAPPING_FILE="00.samples.final.txt"

# 最终输出的汇总表格
OUTPUT_TABLE="final_summary_table.tsv"

# 用于分析的CPU线程数
THREADS=8

# BUSCO谱系库的路径 (请根据您的BUSCO安装和物种类型修改)
# 例如: embryophyta_odb10 (陆生植物), viridiplantae_odb10 (绿色植物)
BUSCO_LINEAGE="embryophyta_odb10"

# 所有结果的输出总目录
RESULTS_DIR="analysis_results"

# --- SCRIPT SETUP ---

# 如果任何命令失败，脚本将立即退出
set -e
# 如果管道中的任何一个命令失败，则整个管道视为失败
set -o pipefail

# --- DEPENDENCY CHECK ---
# 检查所需软件是否都存在
echo "--> Checking for required software..."
check_dependency() {
    if ! command -v "$1" &> /dev/null; then
        echo "Error: Required software '$1' not found. Please install it or activate the correct conda environment."
        exit 1
    fi
}
check_dependency bwa
check_dependency samtools
check_dependency seqtk
check_dependency busco
echo "All software found."

# --- INITIALIZE ---

# 创建总的结果输出目录
mkdir -p ${RESULTS_DIR}

# 初始化输出文件，并写入更完整的表头
echo -e "SampleID\tAvg_Coverage_Depth\tGenome_Length_bp\tGC_Content_%\tBUSCO_Complete_%\tNum_PCGs\tNum_tRNAs" > ${OUTPUT_TABLE}


# --- MAIN LOOP ---
# 逐行读取映射文件
while IFS=$'\t' read -r sample_id reads1_path reads2_path assembly_path genbank_path
do
    # 忽略注释行或空行
    if [[ "$sample_id" == \#* ]] || [[ -z "$sample_id" ]]; then
        continue
    fi

    echo "======================================================"
    echo "Processing Sample: ${sample_id}"
    echo "======================================================"

    # 为每个样本创建独立的输出目录
    SAMPLE_DIR="${RESULTS_DIR}/${sample_id}"
    mkdir -p ${SAMPLE_DIR}

    # --- 1. 计算平均测序深度 ---
    echo "--> Calculating Coverage Depth..."
    bwa index "${assembly_path}"
    bwa mem -t ${THREADS} "${assembly_path}" "${reads1_path}" "${reads2_path}" > "${SAMPLE_DIR}/alignment.sam"
    samtools view -bS "${SAMPLE_DIR}/alignment.sam" | samtools sort -@ ${THREADS} -o "${SAMPLE_DIR}/sorted.bam" -
    rm "${SAMPLE_DIR}/alignment.sam" # 删除巨大的sam文件
    coverage=$(samtools depth "${SAMPLE_DIR}/sorted.bam" | awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }')
    printf "Coverage for ${sample_id} is %.2f\n" ${coverage}

    # --- 2. 计算基因组长度和GC含量 ---
    echo "--> Calculating Length and GC Content..."
    seqtk_output=$(seqtk comp "${assembly_path}")
    length=$(echo "$seqtk_output" | awk 'NR>1 {print $2}')
    gc_content=$(echo "$seqtk_output" | awk 'NR>1 {print $5}')
    echo "Length: ${length} bp, GC: ${gc_content}%"
   
    # --- 4. 统计注释基因数量 ---
    echo "--> Calculating Annotation Stats..."
    # 从GenBank文件中用grep统计CDS和tRNA的数量
    num_pcgs=$(grep -c '     CDS     ' "${genbank_path}")
    num_trnas=$(grep -c '     tRNA    ' "${genbank_path}")
    echo "Found ${num_pcgs} PCGs and ${num_trnas} tRNAs."

    # --- 5. 汇总结果并写入输出文件 ---
    echo "--> Writing results to table..."
    echo -e "${sample_id}\t${coverage}\t${length}\t${gc_content}\t${busco_score}\t${num_pcgs}\t${num_trnas}" >> ${OUTPUT_TABLE}

done < ${MAPPING_FILE}

echo "======================================================"
echo "All 62 samples processed successfully!"
echo "Final results table saved to: ${OUTPUT_TABLE}"
echo "======================================================"
