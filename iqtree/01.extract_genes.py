import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord # <-- 新增了这一行导入
from collections import defaultdict

# --- 用户配置 ---
# 存放所有GenBank文件的文件夹路径
GENBANK_DIR = "/data01/wangkun/projects/taibaishaolan/03.analysis/genbank/"
# 输出每个基因FASTA文件的文件夹路径
OUTPUT_DIR = "individual_gene_fastas/"
# -----------------

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 使用defaultdict来自动处理新的基因
genes = defaultdict(list)
sample_files = [f for f in os.listdir(GENBANK_DIR) if f.endswith((".gb", ".gbk"))]

print(f"找到 {len(sample_files)} 个GenBank文件，开始提取基因...")

# 遍历每个GenBank文件
for filename in sample_files:
    filepath = os.path.join(GENBANK_DIR, filename)
    sample_id = os.path.splitext(filename)[0]
    
    record = SeqIO.read(filepath, "genbank")
    
    # 遍历该文件中的所有特征
    for feature in record.features:
        # 我们只关心蛋白质编码基因 (CDS)
        if feature.type == "CDS":
            # 获取基因名，如果没有则跳过
            if "gene" in feature.qualifiers:
                gene_name = feature.qualifiers["gene"][0]
                
                # --- *** 这里是修正的部分 *** ---
                # 步骤 1: 先提取序列对象 (这部分是正确的)
                extracted_seq = feature.extract(record.seq)
                # 步骤 2: 使用正确的 SeqRecord 类来创建一个新的、带名字的序列记录
                seq_record = SeqRecord(extracted_seq, id=sample_id, description="")
                # ---------------------------------
                
                # 将这个序列记录添加到对应的基因字典中
                genes[gene_name].append(seq_record)

print("基因提取完成，开始写入文件...")

# 遍历字典，为每个基因写入一个FASTA文件
written_count = 0
# 使用 sorted() 来确保每次运行输出的文件顺序一致
for gene_name in sorted(genes.keys()):
    records = genes[gene_name]
    # 通常我们只对比对所有（或绝大多数）样本都存在的基因感兴趣
    if len(records) >= len(sample_files) * 0.9: # 例如，只保留在90%以上样本中都存在的基因
        output_filename = os.path.join(OUTPUT_DIR, f"{gene_name}.fasta")
        SeqIO.write(records, output_filename, "fasta")
        written_count += 1

print(f"\n处理完成！共为 {written_count} 个共享基因创建了FASTA文件。")
print(f"文件保存在: {OUTPUT_DIR}")
