import os
import pandas as pd
from time import time
from Bio import SeqIO
import math

def positional_encode(sequence, alpha=0.05):

    mapping = {'A':0, 'C':1, 'G':2, 'T':3}
    encoded_matrix = []

    length = len(sequence)
    count_A = sequence.count('A')
    count_C = sequence.count('C')
    count_G = sequence.count('G')
    count_T = sequence.count('T')

    cls_token = [
        count_A / length,
        count_C / length,
        count_G / length,
        count_T / length
    ]

    encoded_matrix.append(cls_token)
    
    for k_1based, nt in enumerate(sequence, 1):
        if nt not in mapping:
            continue

        vec = [0.0] * 4
        vec[mapping[nt]] = 1.0

        angle = 2 * math.pi * k_1based / 10.0
        sin_val = math.sin(angle)
        cos_val = math.cos(angle)

        pos_enc = [sin_val, cos_val, -sin_val, cos_val]

        fused = [v + alpha*p for v,p in zip(vec, pos_enc)]
        encoded_matrix.append(fused)

    return encoded_matrix


def process_file(input_file, output_base_dir):
  
    file_name = os.path.splitext(os.path.basename(input_file))[0]
    
    output_folder = os.path.join(output_base_dir, file_name)
    
    os.makedirs(output_folder, exist_ok=True)
    
    sequences = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        seq = str(record.seq).upper()

    columns = ['A', 'C', 'G', 'T']
    
    count = 0
    for name, seq in sequences:
        encoded = positional_encode(seq)
        if encoded:
            df = pd.DataFrame(encoded, columns=columns)
            safe_name = "".join(c for c in name if c.isalnum() or c in ('_', '-')).rstrip()
            output_file = os.path.join(output_folder, f"{safe_name}.csv")
            df.to_csv(output_file, index=False)
            count += 1
    
    return count, len(sequences)


def process_all_txt_files(input_folder, output_base_dir):
    start_time = time()
    
    os.makedirs(output_base_dir, exist_ok=True)
    
    txt_files = [f for f in os.listdir(input_folder) if f.endswith('.txt')]
    
    if not txt_files:
        print("not found")
        return

    for i, txt_file in enumerate(txt_files, 1):
        input_file_path = os.path.join(input_folder, txt_file)
        
        seq_processed, seq_found = process_file(
            input_file_path, output_base_dir
        )
            
    end_time = time()
    

def main():
    # 配置路径
    input_folder = "/data/test"  # 包含所有txt文件的文件夹
    output_base_dir = "/data/test/mfpe"  # 输出基础目录
    
    if not os.path.exists(input_folder):
        print(f"error: {input_folder}")
        return
    
    process_all_txt_files(input_folder, output_base_dir)
    



if __name__ == "__main__":
    main()
