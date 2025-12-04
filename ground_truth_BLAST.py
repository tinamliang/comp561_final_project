from collections import defaultdict
import pandas as pd

"""
find overlapping reads based on their read-to-genome positions from a BLAST output file. 
where suffix of Ri matches with prefix from Rj.
"""

def get_start_end_pos(blast_f):
    results = defaultdict(list)
    with open(blast_f, 'r') as f:
        for line in f:
            line = line.strip()
            if line :
                fields = line.split('\t')
                results[fields[0]].append((int(fields[1]), int(fields[2])))
    return results

def intervals_overlap(int1, int2):
    start_idx = max(int1[0], int2[0])
    end_idx = min(int1[1], int2[1])

    if end_idx > start_idx: 
        return (end_idx - start_idx)
    else:
        return 0

def reads_overlap(intervals1, intervals2, read1_id, read2_id):
    """Check if any interval from read1 overlaps with any interval from read2"""
    overlaps = []
    for i1 in intervals1:
        for i2 in intervals2:
            interval = intervals_overlap(i1, i2)
            if interval > 100:
                overlaps.append((read1_id, read2_id, i1, i2))
    return overlaps

def find_overlapping_reads(reads):
    """Find all pairs of overlapping reads"""
    overlap_list = []
    read_ids = list(reads.keys())
    
    for i in range(len(read_ids)):
        for j in range(i + 1, len(read_ids)):
            ri, rj = read_ids[i], read_ids[j]
            overlaps = reads_overlap(reads[ri], reads[rj], ri, rj)
            overlap_list.extend(overlaps)
    print(len(overlap_list))
    return overlap_list

def overlaps_to_table(overlaps, output_file='overlaps_ground_truth.csv'):
    """Convert overlap tuples to a table with chrom_id as left column"""
    rows = []
    for overlap in overlaps:
        read1, read2, idx1, idx2 = overlap
        rows.append({
            'read1_id': read1,
            'read2_id': read2,
            'read1_start': idx1[0],
            'read1_end': idx1[1],
            'read2_start': idx2[0],
            'read2_end': idx2[1],
        })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Table written to {output_file}")
    return df

if __name__ == "__main__":
    # this file has reads aligning to the area in the genome with custom output format
    blast_file = 'alignments_custom_tab.txt' 
    results = get_start_end_pos(blast_file)
    overlaps = find_overlapping_reads(results)
    df = overlaps_to_table(overlaps)