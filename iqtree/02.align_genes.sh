# 创建一个存放比对结果的新目录
mkdir -p ../aligned_genes

# 循环处理当前目录下的每一个.fasta文件
for file in *.fasta
do
    echo "Aligning ${file}..."
    # 使用MAFFT进行比对，并将结果保存到新目录中，后缀为.aln.fasta
    mafft --auto --thread 8 "${file}" > "../aligned_genes/${file%.fasta}.aln.fasta"
done
