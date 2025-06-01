from Bio import SeqIO # 使用SeqIO模块来处理fastq文件
import numpy as np
from collections import Counter
import time
import argparse
import os
import tarfile


def output_non_sequence_info(id: list, length: list, quality: list, output_filepath: str):
    """
    输出非序列信息，包括ID、长度和质量。
    """
    try:
        with open(output_filepath, 'w', encoding='utf-8') as outfile:
            leng = len(id)
            for i in range(leng):
                outfile.write(f"@{id[i]}\n{length[i]}\n{quality[i]}\n")
    except Exception as e:
        print(f"Error when writing '{output_filepath}': {e}")


def _count_overlapping(text: str, pattern: str) -> int:
    """
    高效统计子字符串 pattern 在 text 中重叠出现的次数。
    """
    count = 0
    i = -1 # 从文本的起始位置开始查找
    text_len = len(text)
    pattern_len = len(pattern)

    if pattern_len == 0: # 空字符串的特殊处理
        return text_len + 1 # 通常定义为空字符串在文本中出现 len(text) + 1 次

    if pattern_len > text_len: # 模式串比文本长，不可能出现
        return 0

    while True:
        i = text.find(pattern, i + 1) # 从上一个找到的位置之后继续查找
        if i == -1: # 未找到更多匹配
            break
        count += 1
    return count


def get_kmer(k: int, seq: str):
    """
    按照原序列的碱基比例随机生成k-mer，然后找到这些生成的k-mer中，在原序列里重复次数最高的那个。
    """
    # 根据当前时间设置唯一的随机种子
    np.random.seed(int(time.time())) # 使用 time.time() 获取浮点数时间戳

    len_seq = len(seq)

    if k <= 0:
        raise ValueError("k must be positive.")

    if len_seq < k: # 此条件仅当 k > 0 时有意义
        raise ValueError("Sequence length must be at least k.")
    
    # Step 1: 计算碱基频率
    # 仅考虑 'A', 'C', 'G', 'T', 'N'
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
    for base in seq:
        if base in counts:
            counts[base] += 1

    proportions = {base: count / len_seq for base, count in counts.items()}
    
    # Step 2: 根据碱基比例生成100个随机k-mer
    bases_for_choice = list(proportions.keys())
    probs_for_choice = np.array(list(proportions.values()))
    
    # 确保概率总和为1 (处理浮点数精度问题)
    prob_sum = probs_for_choice.sum()
    probs_for_choice /= prob_sum # 归一化
    
    kmers_generated = []
    for _ in range(100): # 生成kmers
        kmer = ''.join(np.random.choice(bases_for_choice, size=k, p=probs_for_choice))
        kmers_generated.append(kmer)
    
    # Step 3: 统计每个生成的k-mer在原始序列中的出现次数
    kmer_counts = Counter()
    unique_generated_kmers = set(kmers_generated) # 只对不重复的生成k-mer进行计数

    for kmer_to_count in unique_generated_kmers:
        # 使用优化后的计数函数
        count = _count_overlapping(seq, kmer_to_count)
        kmer_counts[kmer_to_count] = count

    # 获取出现次数最多的k-mer及其次数
    # kmer_counts.most_common(1) 返回一个列表，如 [('ACG', 5)]
    # 如果 kmer_counts 为空，这里会抛出 IndexError
    most_common_kmer, count = kmer_counts.most_common(1)[0]
    
    return most_common_kmer


# 下面为针对已经完成kmer压缩后的类二进制序列压缩算法

# --- 碱基对与字符的映射定义 ---
# 这个映射在压缩和解压缩时必须完全一致

BASES = ['A', 'C', 'G', 'T', 'N']
# 我们需要 5*5 = 25 个唯一字符来映射所有碱基对。
# 这里我们使用小写字母 'a' 到 'y'。

MAPPING_CHARS = [chr(ord('a') + i) for i in range(25)] # 'a' 到 'y'

DNA_PAIR_TO_CHAR = {}
CHAR_TO_DNA_PAIR = {}

_char_iterator = iter(MAPPING_CHARS)
for base1 in BASES:
    for base2 in BASES:
        pair = base1 + base2
        char_representation = next(_char_iterator)
        DNA_PAIR_TO_CHAR[pair] = char_representation
        CHAR_TO_DNA_PAIR[char_representation] = pair
# --- 定义结束 ---


def compress_dna_to_file(input_filepath: str, output_filepath: str):
    """
    压缩DNA序列文件。
    读取DNA序列，将每两个碱基替换为一个映射字符。
    如果序列长度为奇数，最后一个碱基保持原样。
    """
    try:
        with open(input_filepath, 'r', encoding='utf-8') as infile:
            dna_sequence = infile.readline().strip() # 读取并去除末尾换行符
    except FileNotFoundError:
        print(f"Error: Input file '{input_filepath}' not found.")
        return
    except Exception as e:
        print(f"Error in reading '{input_filepath}' : {e}")
        return

    # 验证输入序列是否只包含允许的碱基
    for base in dna_sequence:
        if base not in BASES:
            print(f"Error: Invalid base '{base}' in DNA Sequence。Allowed bases are {BASES}.")
            return

    compressed_parts = []
    i = 0
    n = len(dna_sequence)
    while i < n:
        if i + 1 < n:  # 如果可以凑成一对
            pair = dna_sequence[i:i+2]
            try:
                compressed_parts.append(DNA_PAIR_TO_CHAR[pair])
            except KeyError:
                # 如果DNA序列中包含非ACGTN字符导致无法形成有效pair
                print(f"Error: Invalid DNA Pair '{pair}'")
                exit(1)

            i += 2
        else:  # 序列长度为奇数，剩下最后一个碱基
            compressed_parts.append(dna_sequence[i])
            i += 1
            
    compressed_string = "".join(compressed_parts)

    try:
        with open(output_filepath, 'w', encoding='utf-8') as outfile:
            outfile.write(compressed_string)
    except Exception as e:
        print(f"Error occurrence when writing '{output_filepath}': {e}")


def decompress_dna_from_file(input_filepath: str, output_filepath: str):
    """
    解压缩由 compress_dna_to_file 函数生成的压缩文件，还原原始DNA序列。
    """
    try:
        with open(input_filepath, 'r', encoding='utf-8') as infile:
            compressed_sequence = infile.readline().strip()
    except FileNotFoundError:
        print(f"Error: Input file '{input_filepath}' not found.")
        return
    except Exception as e:
        print(f"Error when reading '{input_filepath}': {e}")
        return

    decompressed_parts = []
    for char_code in compressed_sequence:
        if char_code in CHAR_TO_DNA_PAIR: # 这是一个被压缩的碱基对
            decompressed_parts.append(CHAR_TO_DNA_PAIR[char_code])
        elif char_code in BASES: # 这是原始序列中未被压缩的单个末尾碱基
            decompressed_parts.append(char_code)
        else:
            # 如果字符既不是映射字符，也不是有效的单个碱基
            # 这可能表示压缩文件已损坏或包含无效字符
            print(f"WARNING: Unrecognizable character '{char_code}' in compressed file. It will be kept as is.")
            decompressed_parts.append(char_code) # 原样保留
    decompressed_string = "".join(decompressed_parts)

    try:
        with open(output_filepath, 'w', encoding='utf-8') as outfile:
            outfile.write(decompressed_string)
    except Exception as e:
        print(f"Error when writing '{output_filepath}': {e}")


def cut_sequence_by_length(seq: str, length: list) -> list:
    """
    将序列按指定长度切分成多个片段。
    如果最后一个片段不足指定长度，则保留原样。
    """
    n = len(length)
    result = []
    pos = 0
    for i in range(n):
        result.append(seq[pos:pos + length[i]])
        pos += length[i]
    return result


def kmerdeletion(seq,most_frequent_kmer):
    positions=[]
    k=len(most_frequent_kmer)
    i=0
    
    # 记录重复序列出现的位置
    while i<=len(seq)-k:
        if seq[i:i+k]==most_frequent_kmer:
            positions.append(i)
            i+=k
        else:
            i+=1
    
    # 记录非重复序列的位置      
    segments = []
    last_end = 0
    for pos in positions:
        segments.append(seq[last_end:pos])  # 添加k-mer前的片段
        last_end = pos + k  # 跳过k-mer
    
    segments.append(seq[last_end:])  # 添加最后一段
    new_seq = "".join(segments)
    
    # 利用差值保存重复位点信息
    positions=[positions[0]]+[positions[i]-positions[i-1] for i in range(1,len(positions))]
    
    return new_seq,positions


def kmerrecover(seq_file, kmer_pos_file):
    # 读取删除后的序列
    with open(seq_file,"r") as f:
        new_seq=""
        for line in f.readlines():
            line.strip()
            new_seq+=line
    
    # 读取kmer和位置信息
    with open(kmer_pos_file) as f:
        most_frequent_kmer=""
        positions=[]
        for line in f.readlines():
            line=line.strip()
            if line and line[0].isalpha():
                most_frequent_kmer+=line
            else:
                positions=line.split(",")
                positions.pop()
                positions=list(map(int,positions))

          
    # 将重复位置信息复原
    positions=[sum(positions[:i+1]) for i in range(len(positions))]
    
    # 复原序列
    for pos in positions:
        if pos<len(new_seq):
            new_seq=new_seq[:pos]+most_frequent_kmer+new_seq[pos:]
        else:
            new_seq+=most_frequent_kmer
    return new_seq


def seq_output(output_file, new_seq):
    with open(output_file,"w") as f:
        f.write(f"{new_seq}\n")


def kmer_pos_output(output_file, kmer, positions):
    with open(output_file,"w") as f:
        f.write(f"{kmer}\n")
        for pos in positions:
            f.write(f"{pos},")


def hamming_distance(seq1, seq2):
    return sum(b1 != b2 for b1, b2 in zip(seq1, seq2))


def find_all_matches(sequence, seed, max_mismatches=1): 
    """ 查找所有符合条件的子序列，返回错配位置索引 """
    seed_length = len(seed)  
    matches = []
    for i in range(len(sequence) - seed_length + 1):
        subsequence = sequence[i:i+seed_length]
        distance = hamming_distance(subsequence, seed)
        if distance <= max_mismatches:
            mismatch_positions = [pos for pos, (b1, b2) in enumerate(zip(subsequence, seed)) if b1 != b2]
            matches.append({
                'start': i,
                'subsequence': subsequence,
                'hamming_distance': distance,
                'mismatch_positions': mismatch_positions
            })
    return matches


def kmerdeletion_optimized(file_path, seed, output_base_path, max_mismatches=1):  
    seed_length = len(seed)
    with open(file_path, 'r') as f:
        sequence = f.read().strip()
        clean_sequence = ''.join([base for base in sequence if base in 'ATCGN'])
        if len(clean_sequence) < seed_length:
            print(f"kmer length {seed_length} is longer than sequence.")
            return
        
        all_matches = find_all_matches(clean_sequence, seed, max_mismatches)
        if not all_matches:
            print(f"Unable to find sub sequences with mismatches≤{max_mismatches}")
            return
        
        # 正序排序以计算正确的相对位置
        all_matches_sorted = sorted(all_matches, key=lambda x: x['start'])
        
        deleted_info = []
        prev_end = 0
        for match in all_matches_sorted:
            start = match['start']
            end = start + seed_length
            relative_start = start - prev_end
            
            mismatch_positions = match['mismatch_positions']
            if not mismatch_positions:
                formatted_mismatch = 'n'
                formatted_bases = 'n'
            else:
                formatted_mismatch = ';'.join(map(str, mismatch_positions))
                formatted_bases = ';'.join([match['subsequence'][pos] for pos in mismatch_positions])
            
            deleted_info.append({
                'relative_start': relative_start,
                'mismatch_pos': formatted_mismatch,
                'original_bases': formatted_bases,
                'absolute_start': start  # 仅用于删除，不保存到文件
            })
            prev_end = end

        # 实际删除子序列：倒序删除避免位置错乱
        clean_list = list(clean_sequence)
        for match in sorted(deleted_info, key=lambda x: x['absolute_start'], reverse=True):
            del clean_list[match['absolute_start']:match['absolute_start'] + seed_length]

        processed_sequence = ''.join(clean_list)
        total_deleted_bases = len(clean_sequence) - len(processed_sequence)

        # 写入处理后的序列
        with open(f"{output_base_path}.seqtemp", 'w') as f:
            f.write(processed_sequence)
        
        # 写入删除信息文件
        with open(f"{output_base_path}.pos", 'w') as f:
            f.write(f"{seed}\n")
            for item in deleted_info:
                f.write(f"{item['relative_start']},{item['mismatch_pos']},{item['original_bases']}\n")


def kmerrecover_optimized(processed_file_path, delete_info_file_path):
    with open(processed_file_path, 'r') as f:
        processed_seq = list(f.read().strip())

    with open(delete_info_file_path, 'r') as f:
        seed = f.readline().strip()
        seed_len = len(seed)
        deletion_info = [line.strip().split(',') for line in f if line.strip()]

    prev_end = None  # 用于计算相对偏移，None 表示第一个是绝对位置

    for rel_start_str, mismatch_pos, original_bases in deletion_info:
        rel_start = int(rel_start_str)

        # 计算插入位置
        if prev_end is None:
            absolute_start = rel_start  # 第一个是绝对位置
        else:
            absolute_start = prev_end + rel_start  # 其余是相对偏移

        # 构建恢复的子序列
        subsequence = list(seed)
        if mismatch_pos != 'n':
            pos_list = list(map(int, mismatch_pos.split(';')))
            base_list = original_bases.split(';') if original_bases != 'n' else []
            for pos, base in zip(pos_list, base_list):
                subsequence[pos] = base
        subsequence = ''.join(subsequence)

        # 插入子序列
        processed_seq[absolute_start:absolute_start] = subsequence

        # 更新 prev_end
        prev_end = absolute_start + seed_len
    
    return processed_seq

def extract_targz(targz_file):
    """
    Extract a TAR.GZ file to the specified directory.
    
    Args:
        targz_file (str): Path to the input TAR.GZ file
        output_dir (str): Directory to extract files to (default: current directory)
    """
    try:
        # Check if the file exists
        if not os.path.isfile(targz_file):
            print(f"Error: {targz_file} does not exist or is not a file")
            return
        output_dir = os.path.dirname(os.path.abspath(targz_file))

        if not output_dir:  # 如果文件在当前目录，dirname 返回空字符串
            output_dir = os.getcwd()
        
        # Open and extract the TAR.GZ file
        with tarfile.open(targz_file, 'r:gz') as tarf:
            tarf.extractall(path=output_dir, filter='data')
                
    except tarfile.TarError as e:
        print(f"Error extracting {targz_file}: {e}")
    except OSError as e:
        print(f"Error accessing files: {e}")


def parse_file_to_lists(input_file):
    """
    Parse an input file where each group of three lines corresponds to ID, length, and quality.
    Store the data in three lists: id, length, and quality.

    Args:
        input_file (str): Path to the input file

    Returns:
        tuple: Three lists (id, length, quality)
    """
    id_list = []
    length_list = []
    quality_list = []

    try:
        with open(input_file, 'r') as file:
            lines = file.readlines()
            # Remove any trailing newlines
            lines = [line.strip() for line in lines]

            # Check if the number of lines is a multiple of 3
            if len(lines) % 3 != 0:
                print(f"Warning: File {input_file} has {len(lines)} lines, which is not a multiple of 3")
                return id_list, length_list, quality_list

            # Process lines in groups of 3
            for i in range(0, len(lines), 3):
                id_list.append(lines[i])
                try:
                    length_list.append(int(lines[i + 1]))  # Convert length to integer
                except ValueError:
                    print(f"Error: Length at line {i + 2} ('{lines[i + 1]}') is not an integer")
                    return id_list, length_list, quality_list
                quality_list.append(lines[i + 2])

        return id_list, length_list, quality_list

    except FileNotFoundError:
        print(f"Error: File {input_file} not found")
    except Exception as e:
        print(f"Error reading file {input_file}: {e}")


def write_fastq(output_file, id_list, sequence_list, quality_list):
    """
    Write ID, sequence, and quality lists to a FASTQ file.

    Args:
        output_file (str): Path to the output FASTQ file
        id_list (list): List of sequence IDs
        sequence_list (list): List of sequences
        quality_list (list): List of quality score strings
    """
    # Check if lists have the same length
    if not (len(id_list) == len(sequence_list) == len(quality_list)):
        print(f"Error: Lists have different lengths: "
              f"IDs ({len(id_list)}), Sequences ({len(sequence_list)}), "
              f"Qualities ({len(quality_list)})")
        return

    try:
        with open(output_file, 'w') as f:
            for idx, (id_, seq, qual) in enumerate(zip(id_list, sequence_list, quality_list)):
                # Write FASTQ record
                f.write(f"{id_}\n")
                f.write(f"{seq}\n")
                f.write(f"+\n")
                f.write(f"{qual}\n")

    except OSError as e:
        print(f"Error writing to {output_file}: {e}")


def compress_strategy1(inputfile, outputfile):
    """
    压缩策略1：使用k-mer压缩。
    读取输入文件，找到最频繁的k-mer，并删除其所有出现位置。
    输出压缩后的序列和k-mer位置。
    """
    
    fastq_file = inputfile
    id = []
    length = []
    quality = []
    sequences = [] # 初始化四个列表来存储fastq文件的内容

    print(f"Reading FASTQ file: {fastq_file}")
    for record in SeqIO.parse(fastq_file, "fastq"):
        fastq_lines = record.format("fastq").strip().split("\n")
        id.append(record.id) # 记录ID
        length.append(len(record.seq)) # 记录序列长度
        sequences.append(str(record.seq)) # 记录序列
        quality_string = fastq_lines[3] # 记录质量字符串
        quality.append(quality_string)

    combined_sequences = "".join(sequences)

    for base in combined_sequences: # 检查序列是否合规
        if base not in "ACGTN":
            raise ValueError(f"Found Invalid Character in sequences: {base}")

    output_non_sequence_info(id, length, quality, outputfile + ".info")  # 输出非序列信息
    print("Non-sequence information has been saved. ")
    print("Starting standard k-mer compression...")
    kmer = get_kmer(6, combined_sequences) # 获得最频繁的kmer
    new_seq, positions = kmerdeletion(combined_sequences, kmer) # 实行kmer删除
    print("Standard k-mer compression completed.")
    # 输出压缩后的序列和k-mer位置
    seq_output(outputfile + ".seqtemp", new_seq)
    kmer_pos_output(outputfile + ".pos", kmer, positions)
    print("kmer positions have been saved.")
    # 进行类二进制压缩序列文件
    compress_dna_to_file(outputfile + ".seqtemp", outputfile + ".seq")  # 压缩序列文件
    os.remove(outputfile + ".seqtemp")
    print("Compressed sequence has been saved. Now starting packaging.")
    try: # 用tar-gzip进一步打包，压缩文件
        with tarfile.open(outputfile + ".tar.gz", 'w:gz') as tarf:
            files = [outputfile + ".pos", outputfile + ".seq", outputfile + ".info"]
            for file_path in files:
                if os.path.isfile(file_path):
                    tarf.add(file_path, arcname=os.path.basename(file_path))
                    print(f"Added: {file_path}")
                else:
                    print(f"Skipped: {file_path} does not exist or is not a file")
        print(f"Created TAR.GZ archive: {outputfile + ".tar.gz"}")
        os.remove(outputfile + ".seq")
        os.remove(outputfile + ".pos")
        os.remove(outputfile + ".info")
    except Exception as e:
        print(f"Error creating TAR.GZ archive: {e}")
    print("Compression strategy 1 completed successfully.")


def compress_strategy2(inputfile, outputfile):
    """
    压缩策略2：使用优化的k-mer压缩。
    读取输入文件，找到最频繁的k-mer，并删除其所有出现位置与相似序列出现位置。
    输出压缩后的序列和k-mer位置。
    """
    
    fastq_file = inputfile
    id = []
    length = []
    quality = []
    sequences = [] # 初始化四个列表来存储fastq文件的内容

    print(f"Reading FASTQ file: {fastq_file}")
    for record in SeqIO.parse(fastq_file, "fastq"):
        fastq_lines = record.format("fastq").strip().split("\n")
        id.append(record.id) # 记录ID
        length.append(len(record.seq)) # 记录序列长度
        sequences.append(str(record.seq)) # 记录序列
        quality_string = fastq_lines[3] # 记录质量字符串
        quality.append(quality_string)

    combined_sequences = "".join(sequences)
    
    for base in combined_sequences: # 检查序列是否合规
        if base not in "ACGTN":
            raise ValueError(f"Found Invalid Character in sequences: {base}")
    
    with open(f"{outputfile}.temp", 'w') as f:
        f.write(combined_sequences)

    output_non_sequence_info(id, length, quality, outputfile + ".info")  # 输出非序列信息
    print("Non-sequence information has been saved. ")
    print("Starting optimized k-mer compression...")
    kmer = get_kmer(6, combined_sequences) # 获得最频繁的kmer
    kmerdeletion_optimized(outputfile + ".temp", kmer, outputfile) # 实行kmer删除
    print("Optimized k-mer compression completed.")
    print("kmer positions have been saved.")
    # 进行类二进制压缩序列文件
    compress_dna_to_file(outputfile + ".seqtemp", outputfile + ".seq")  # 压缩序列文件
    os.remove(outputfile + ".seqtemp")
    print("Compressed sequence has been saved. Now starting packaging.")
    try: # 用tar-gzip进一步打包，压缩文件
        with tarfile.open(outputfile + ".tar.gz", 'w:gz') as tarf:
            files = [outputfile + ".pos", outputfile + ".seq", outputfile + ".info"]
            for file_path in files:
                if os.path.isfile(file_path):
                    tarf.add(file_path, arcname=os.path.basename(file_path))
                    print(f"Added: {file_path}")
                else:
                    print(f"Skipped: {file_path} does not exist or is not a file")
        print(f"Created TAR.GZ archive: {outputfile + ".tar.gz"}")
        os.remove(outputfile + ".seq")
        os.remove(outputfile + ".pos")
        os.remove(outputfile + ".info")
        os.remove(outputfile + ".temp")
    except Exception as e:
        print(f"Error creating TAR.GZ archive: {e}")
    print("Compression strategy 2 completed successfully.")


def decompress_strategy1(inputfile, outputfile):
    """
    将使用压缩策略1的文件解压
    """
    basename = os.path.splitext(inputfile)[0]
    basename = os.path.splitext(basename)[0]  # 移除 .tar.gz后缀名
    print("Extracting the input file.")
    extract_targz(inputfile) # 解压文件到inpufile目录（一会会删除）
    inputfile = basename
    print("Decompressing the sequence file.")
    decompress_dna_from_file(inputfile + ".seq", inputfile + ".seqtemp")
    print("Decompressing the kmer. ")
    tempseq = kmerrecover(inputfile + ".seqtemp", inputfile + ".pos")
    print("Getting the other information.")
    id, length, quality = parse_file_to_lists(inputfile + ".info")
    sequence = cut_sequence_by_length(tempseq, length)
    print("Writing the output file.")
    write_fastq(outputfile, id, sequence, quality)
    os.remove(inputfile + ".seq")
    os.remove(inputfile + ".seqtemp")
    os.remove(inputfile + ".pos")
    os.remove(inputfile + ".info")
    print("Decompression Method 1 Completed")
    
    
def decompress_strategy2(inputfile, outputfile):
    """
    将使用压缩策略2的文件解压
    """
    basename = os.path.splitext(inputfile)[0]
    basename = os.path.splitext(basename)[0]  # 移除 .tar.gz后缀名
    print("Extracting the input file.")
    extract_targz(inputfile) # 解压文件到inpufile目录（一会会删除）
    inputfile = basename
    print("Decompressing the sequence file.")
    decompress_dna_from_file(inputfile + ".seq", inputfile + ".seqtemp")
    print("Decompressing the kmer. ")
    tempseq = "".join(kmerrecover_optimized(inputfile + ".seqtemp", inputfile + ".pos"))
    print("Getting the other information.")
    id, length, quality = parse_file_to_lists(inputfile + ".info")
    sequence = cut_sequence_by_length(tempseq, length)
    print("Writing the output file.")
    write_fastq(outputfile, id, sequence, quality)
    os.remove(inputfile + ".seq")
    os.remove(inputfile + ".seqtemp")
    os.remove(inputfile + ".pos")
    os.remove(inputfile + ".info")
    print("Decompression Method 2 Completed")


def main():
    parser = argparse.ArgumentParser(description="A Simple FASTQ Compressor and Decompressor")

    parser.add_argument('--type', type=int, required=True, choices=[1, 2],
                        help='1 = compress, 2 = decompress')
    parser.add_argument('--method', type=int, choices=[1, 2],
                        help='Compression method (1-standard kmer compression + binary-like compression, 2-optimized kmer compression + binary-like compression)')
    parser.add_argument('inputfile', type=str, help='Input FASTQ file/ TAR.GZ file')
    parser.add_argument('outputfile', type=str, help='Output file (if compressing, this will be the prefix of a TAR.GZ file; if decompressing, this will be a FASTQ file)')

    args = parser.parse_args()

    if args.type == 1:
        if args.method is None:
            print("--method is not given. Defaulting to method 1.")
            args.method = 1
        if args.method == 1:
            compress_strategy1(args.inputfile, args.outputfile)
        elif args.method == 2:
            compress_strategy2(args.inputfile, args.outputfile)
        else:
            raise ValueError("Unknown compression method")
    elif args.type == 2:
        if args.method is None:
            print("--method is not given. Automatically trying compression method 1.")
            try:
                decompress_strategy1(args.inputfile, args.outputfile)
            except Exception as e:
                print(f"Method 1 failed: {e}. Trying Method 2.")
                decompress_strategy2(args.inputfile, args.outputfile)
        elif args.method == 1:
            decompress_strategy1(args.inputfile, args.outputfile)
        elif args.method == 2:
            decompress_strategy2(args.inputfile, args.outputfile)
        else:
            raise ValueError("Unknown compression method")

if __name__ == '__main__':
    main()
