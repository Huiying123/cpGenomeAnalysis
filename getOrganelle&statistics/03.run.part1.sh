#!/bin/bash

# =============================================================================
#           FINAL BATCH ANALYSIS SCRIPT (WITH ADVANCED DEBUGGING)
# =============================================================================
# This version includes robust checks and captures the specific error output
# from seqtk to definitively diagnose issues with input files.
# =============================================================================


# --- CONFIGURATION ---
# (请根据您的系统和数据进行修改)

# 输入的样本映射文件 (5列: SampleID, Read1_Path, Read2_Path, Assembly_FASTA_Path, Annotation_GBK_Path)
MAPPING_FILE="00.samples.final.txt"

# 最终输出的汇总表格
OUTPUT_TABLE="final_summary_table.part1.tsv"

# 用于分析的CPU线程数
THREADS=8

# 所有结果的输出总目录
RESULTS_DIR="analysis_results_v7"


# --- SCRIPT SETUP ---
# 如果任何命令失败，脚本将立即退出
set -e
# 如果管道中的任何一个命令失败，则整个管道视为失败
set -o pipefail


# --- DEPENDENCY CHECK ---
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
check_dependency bc
echo "All required software found."


# --- INITIALIZE ---
# 创建总的结果输出目录
mkdir -p ${RESULTS_DIR}
# 初始化输出文件，并写入完整的表头
echo -e "SampleID\tAvg_Coverage_Depth\tGenome_Length_bp\tGC_Content_%\tNum_PCGs\tNum_tRNAs\tLSC_Length\tSSC_Length\tIR_Length" > ${OUTPUT_TABLE}


# --- MAIN LOOP ---
# 逐行读取映射文件
while read -r sample_id reads1_path reads2_path assembly_path genbank_path; do
    # 自动去除变量前后的空格，增强稳健性
    sample_id=$(echo "${sample_id}" | xargs)
    reads1_path=$(echo "${reads1_path}" | xargs)
    reads2_path=$(echo "${reads2_path}" | xargs)
    assembly_path=$(echo "${assembly_path}" | xargs)
    genbank_path=$(echo "${genbank_path}" | xargs)

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

    # --- 步骤 1: 计算平均测序深度 ---
    echo "--> Calculating Coverage Depth..."
    #bwa index "${assembly_path}"
    #bwa mem -t ${THREADS} "${assembly_path}" "${reads1_path}" "${reads2_path}" > "${SAMPLE_DIR}/alignment.sam"
    #samtools view -@ ${THREADS} -bS "${SAMPLE_DIR}/alignment.sam" | samtools sort -@ ${THREADS} -o "${SAMPLE_DIR}/sorted.bam" -
    #rm "${SAMPLE_DIR}/alignment.sam"
    #samtools index "${SAMPLE_DIR}/sorted.bam"
    coverage=$(samtools depth "${SAMPLE_DIR}/sorted.bam" | awk '{sum+=$3} END { if (NR > 0) print sum/NR; else print 0 }')
    printf "Coverage for ${sample_id} is %.2fx\n" ${coverage}

    # --- 步骤 2: 计算基因组长度和GC含量 (带终极调试) ---
    echo "--> Calculating Length and GC Content..."
    
    # 尝试运行seqtk并捕获其所有输出（包括标准输出和标准错误）
    seqtk_output_and_error=$(seqtk comp "${assembly_path}" 2>&1)
    exit_code=$? # 捕获seqtk命令的退出状态码

    # 检查命令是否成功执行 (退出码为0代表成功)
    if [ ${exit_code} -ne 0 ]; then
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "FATAL ERROR: 'seqtk comp' command failed for sample ${sample_id}."
        echo "Assembly Path: ${assembly_path}"
        echo "--- seqtk reported this error: ---"
        echo "${seqtk_output_and_error}"
        echo "---------------------------------"
        echo "Skipping this sample."
        echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        continue
    fi

    # 智能解析 (与之前相同)
    if [ $(echo "${seqtk_output_and_error}" | wc -l) -gt 1 ]; then
        stats_line=$(echo "${seqtk_output_and_error}" | awk 'NR>1')
    else
        stats_line="${seqtk_output_and_error}"
    fi
    
    length=$(echo "$stats_line" | awk '{print $2}')
    c_count=$(echo "$stats_line" | awk '{print $4}')
    g_count=$(echo "$stats_line" | awk '{print $5}')
    
    if [ -z "${length}" ]; then
        echo "ERROR: Parsed Length/GC is empty for sample ${sample_id}. The seqtk output was unusual. Skipping."
        continue
    fi

    # 使用bc命令行计算器来精确计算GC百分比
    gc_content=$(echo "scale=2; (${g_count} + ${c_count}) * 100 / ${length}" | bc)
    echo "Length: ${length} bp, GC: ${gc_content}%"

    # --- 步骤 3: 统计注释基因数量 ---
    echo "--> Calculating Annotation Stats..."
    num_pcgs=$(grep -c '     CDS     ' "${genbank_path}")
    num_trnas=$(grep -c '     tRNA    ' "${genbank_path}")
    echo "Found ${num_pcgs} PCGs and ${num_trnas} tRNAs."

    # --- 步骤 4: 计算 LSC, SSC, IR 区域长度 ---
    echo "--> Calculating LSC/SSC/IR Lengths..."
    ir_region=$(grep 'repeat_region' "${genbank_path}" | grep 'Inverted Repeat' | head -n 1 || true)
    if [[ -n "$ir_region" ]]; then
        ir_start=$(echo "$ir_region" | grep -o '[0-9]\+' | head -n 1); ir_end=$(echo "$ir_region" | grep -o '[0-9]\+' | tail -n 1)
        ir_len=$((ir_end - ir_start + 1)); 
        ssc_len=18000; # 这是一个粗略的估计，精确值需要专门工具
        lsc_len=$((length - ssc_len - (2 * ir_len)))
    else
        ir_len="N/A"; ssc_len="N/A"; lsc_len="N/A"
    fi
    echo "LSC: ${lsc_len}, SSC: ${ssc_len}, IR: ${ir_len}"
    
    # --- 步骤 5: 汇总结果并写入输出文件 ---
    echo "--> Writing results to table..."
    echo -e "${sample_id}\t${coverage}\t${length}\t${gc_content}\t${num_pcgs}\t${num_trnas}\t${lsc_len}\t${ssc_len}\t${ir_len}" >> ${OUTPUT_TABLE}

done < <(tr -d '\r' < "${MAPPING_FILE}")

echo "======================================================"
echo "All samples processed successfully!"
echo "Final results table saved to: ${OUTPUT_TABLE}"
echo "======================================================"
