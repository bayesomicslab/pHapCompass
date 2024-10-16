import os
import subprocess

def parse_fragment_file(fragment_file_path):
    fragments = []
    with open(fragment_file_path, 'r') as file:
        for line in file:
            tokens = line.strip().split()
            n_reads = int(tokens[0])
            sequence_id = tokens[1]
            read_data = tokens[2:-1]
            G_sequence = tokens[-1]

            reads = []
            for i in range(0, len(read_data), 2):
                pos = int(read_data[i])
                seq = read_data[i+1]
                reads.append((pos, seq))

            fragment = {
                'n_reads': n_reads,
                'sequence_id': sequence_id,
                'reads': reads,
                'G_sequence': G_sequence
            }
            fragments.append(fragment)

    return fragments

def parse_vcf_file(vcf_file_path):
    vcf_records = []
    with open(vcf_file_path, 'r') as file:
        for line in file: 
            if line.startswith('#'):
                continue

            tokens = line.strip().split('\t')
            scaffold_id = tokens[0]
            pos = int(tokens[1])
            chrom = tokens[2]
            ref = tokens[3]
            alt = tokens[4]
            qual = tokens[5]
            filter_status = tokens[6]
            vcf_record = {
                'scaffold_id': scaffold_id,
                'pos': pos,
                'chrom': chrom,
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'filter_status': filter_status
            }
            vcf_records.append(vcf_record)
    return vcf_records

def merge_reads(reads):
    while len(reads) > 2:
        min_distance = None
        min_pair = None
        for i in range(len(reads)):
            for j in range(i+1, len(reads)):
                pos_i, seq_i = reads[i]
                pos_j, seq_j = reads[j]
                distance = abs(pos_i - pos_j)
                if min_distance is None or distance < min_distance: 
                    min_distance = distance
                    min_pair = (i, j)

        if min_pair is None:
            break  

        i, j = min_pair
        merged_read = merge_two_reads(reads[i], reads[j])
        new_reads = [reads[k] for k in range(len(reads)) if k not in min_pair]
        new_reads.append(merged_read)
        reads = new_reads
    return reads 

def merge_two_reads(read1, read2):
    pos1, seq1 = read1
    pos2, seq2 = read2
    start_pos = min(pos1, pos2)
    end_pos = max(pos1 + len(seq1), pos2 + len(seq2))
    merged_length = end_pos - start_pos
    merged_seq = ['N'] * merged_length 

    for i in range(len(seq1)):
        idx = pos1 - start_pos + i
        merged_seq[idx] = seq1[i]

    for i in range(len(seq2)):
        idx = pos2 - start_pos + i
        if merged_seq[idx] == 'N':
            merged_seq[idx] = seq2[i]
        elif merged_seq[idx] != seq2[i]:
            merged_seq[idx] = 'N'

    merged_seq = ''.join(merged_seq)
    return (start_pos, merged_seq)

def get_max_reference_length(fragments, vcf_records):
    max_length = 0
    for fragment in fragments:
        for pos, seq in fragment['reads']:
            end_pos = pos + len(seq)
            if end_pos > max_length:
                max_length = end_pos
    for variant in vcf_records:
        if variant['pos'] > max_length:
            max_length = variant['pos']
    return max_length

def get_reference_info(fragments, vcf_records):
    max_length = get_max_reference_length(fragments, vcf_records)
    return {'chr1': max_length}

def create_sam_records(fragments, vcf_records):
    sam_records = []
    vcf_dict = {}
    for variant in vcf_records:
        if variant['pos'] not in vcf_dict:
            vcf_dict[variant['pos']] = []
        vcf_dict[variant['pos']].append(variant)
    
    for fragment in fragments:
        sequence_id = fragment['sequence_id']
        reads = fragment['reads']
        reads = merge_reads(reads)
    
        if len(reads) == 1:
            flag = 0 
        elif len(reads) == 2:
            flag = 3
        else: 
            continue
    
        for idx, (pos, seq) in enumerate(reads):
            pos += 1   
    
            if flag == 3:
                if idx == 0:
                    read_flag = 67 
                else:
                    read_flag = 131  
            else: 
                read_flag = flag 
    
            qname = sequence_id
            rname = 'chr1'  
            mapq = 255
            cigar = f'{len(seq)}M'
            rnext = '*'
            pnext = 0
            tlen = 0
            qual = '*'  
    
            overlapping_quals = []
            for i in range(len(seq)):
                variant_pos = pos + i
                if variant_pos in vcf_dict:
                    for variant in vcf_dict[variant_pos]:
                        overlapping_quals.append(variant['qual'])
    
            if overlapping_quals:
                qv_tag = f'QV:Z:{",".join(overlapping_quals)}'
            else:
                qv_tag = ''
    
            sam_record = {
                'QNAME' : qname,
                'FLAG' : read_flag,
                'RNAME' : rname, 
                'POS' : pos,
                'MAPQ' : mapq,
                'CIGAR' : cigar,
                'RNEXT': rnext,
                'PNEXT' : pnext, 
                'TLEN' : tlen,
                'SEQ' : seq,
                'QUAL' : qual,
                'OPTIONAL' : qv_tag
            }
            sam_records.append(sam_record)
    return sam_records

def write_sam_file(sam_file_path, sam_records, ref_info):
    with open(sam_file_path, 'w') as file:
        file.write('@HD\tVN:1.0\tSO:unsorted\n')
        for sn, ln in ref_info.items():
            file.write(f'@SQ\tSN:{sn}\tLN:{ln}\n')
        for record in sam_records:
            fields = [
                record['QNAME'],
                str(record['FLAG']),
                record['RNAME'],
                str(record['POS']),
                str(record['MAPQ']),
                record['CIGAR'],
                record['RNEXT'],
                str(record['PNEXT']),
                str(record['TLEN']),
                record['SEQ'],
                record['QUAL']
            ]
            if 'OPTIONAL' in record and record['OPTIONAL']:
                fields.append(record['OPTIONAL'])
            line = '\t'.join(fields) + '\n'
            file.write(line)

def convert_sam_to_bam(sam_file_path, bam_file_path):
    command = f'samtools view -bS {sam_file_path} -o {bam_file_path}'
    subprocess.run(command, shell=True, check=True)

def sort_sam_file(sam_file_path, sorted_sam_file_path):
    command = f'samtools sort {sam_file_path} -o {sorted_sam_file_path}'
    subprocess.run(command, shell=True, check=True)

def index_bam_file(bam_file_path):
    command = f'samtools index {bam_file_path}'
    subprocess.run(command, shell=True, check=True)

def main(fragment_file_path, vcf_file_path):
    fragments = parse_fragment_file(fragment_file_path)
    vcf_records = parse_vcf_file(vcf_file_path)

    max_ref_length = get_max_reference_length(fragments, vcf_records)

    sam_records = create_sam_records(fragments, vcf_records)

    ref_info = get_reference_info(fragments, vcf_records)

    base_name = os.path.basename(fragment_file_path).replace('.frag.txt', '')

    output_dir = os.path.join(os.path.dirname(vcf_file_path), base_name)

    os.makedirs(output_dir, exist_ok=True)

    sam_file_path = os.path.join(output_dir, f'{base_name}.sam')
    sorted_sam_file_path = os.path.join(output_dir, f'{base_name}_sorted.sam')
    bam_file_path = os.path.join(output_dir, f'{base_name}.bam')
    sorted_bam_file_path = os.path.join(output_dir, f'{base_name}_sorted.bam')
    bam_index_path = f'{sorted_bam_file_path}.bai'

    write_sam_file(sam_file_path, sam_records, ref_info)

    sort_sam_file(sam_file_path, sorted_sam_file_path)

    convert_sam_to_bam(sorted_sam_file_path, sorted_bam_file_path)

    index_bam_file(sorted_bam_file_path)

    print(f'SAM file written to {sorted_sam_file_path}')
    print(f'BAM file written to {sorted_bam_file_path}')
    print(f'BAM index written to {bam_index_path}')


if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print('Usage: python converter.py <fragment_file.frag.txt> <input_file.vcf>')
        sys.exit(1)
    fragment_file_path = sys.argv[1]
    vcf_file_path = sys.argv[2]
    main(fragment_file_path, vcf_file_path)
