import os, re

directory = '/data01/sequencing_data/TaiBaiShaoLan-DNA'
output_dir = '/data01/wangkun/projects/taibaishaolan/02.getorganelle/output3'

commands = []
for filename in os.listdir(directory):
    if filename.endswith('_1.fq.gz'):
        # 找到配对文件名
        file_path = os.path.join(directory, filename)
        paired_filename = filename.replace('_1.fq.gz', '_2.fq.gz')
        paired_file_path = os.path.join(directory, paired_filename)
        
        # 检查配对文件是否存在
        if os.path.exists(paired_file_path):
            print("Found paired files:")
            print(filename.strip(), paired_filename)

            basename = os.path.splitext(os.path.basename(filename))[0]
            match = re.match(r'Unknown_([^_]+)_', basename)
            if match:
                basename = match.group(1)
            subfolder_path = os.path.join(output_dir, basename)
            os.makedirs(subfolder_path, exist_ok=True)

                        
            commands.append(f'cd {subfolder_path}; python /data01/wangkun/software/getorganelle/GetOrganelle-1.7.7.1/get_organelle_from_reads.py -1 {file_path}  -2 {paired_file_path}  -k 21,31,45,65,85,105,127  -t 4 -o output -w 70 -F embplant_pt -s /data01/wangkun/projects/huanghuashaolan/02.getorganelle/ref.fasta > full.log')


            
        else:
            print(f"Error: Paired file {paired_filename} not found for {filename}")

script_name = os.path.basename(__file__)
outfile = f'{script_name}.sh'
with open(outfile, 'w') as output:
    for command in commands:
        output.write(f'{command}\n')
