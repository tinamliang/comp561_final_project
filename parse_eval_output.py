import pandas as pd

def parse_eval_output():
    found_read = 0
    df = pd.read_csv('eval_naive_output.csv')
    df['accuracy'] = df['accuracy'].str.replace('%', '').astype(float)

    a_exists = 0
    for a in df['accuracy']:
        if a != 0.00:
            a_exists += 1

    print(f"Number of non-zero accuracy entries: {a_exists}")
    average = df['accuracy'].mean()
    print(f"Average percentage: {average:.2f}%")

    aligned_reads_sot = parse_blast_sot()
    for read in df['query']:
        if read in aligned_reads_sot and df.loc[df['query'] == read, 'accuracy'].values[0] != 0.00:
            found_read += 1
    print(f"Number of aligned reads found in SOT: {found_read}")

def parse_blast_sot():
    aligned_pairs = []
    df = pd.read_csv('overlaps_ground_truth.csv', sep='\t')
    for read in df['read1_id']:
        if read not in aligned_pairs:
            aligned_pairs.append(read)
    return aligned_pairs

if __name__ == "__main__":
    parse_eval_output()
    #parse_blast_sot()