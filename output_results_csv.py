import os
import time
from collections import defaultdict
import pandas as pd

all_overlapping_pairs = []
def overlapping_pairs_for_read(input_file, output_file):
    query = None
    matches = defaultdict(list)
    match_found = 0
    overlapping_pair = []
    ret_value = []

    with open(input_file, 'r') as f:
        next(f)
        next(f)
        for line in f:
            line = line.strip()

            # Split by any whitespace
            fields = line.split()
            query = fields[0]

            # highly aligned reads (Ri, Rj)
            if fields[1] not in matches[query]:
                matches[query].append(fields[1])

            
    with open(output_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Split by any whitespace
            fields = line.split()
            q = fields[0]

            if query == q:
                # ground truth overlapping reads
                overlapping_pair.append(fields[1])


    for m in matches[query]:
        if m in overlapping_pair:
            # matching Rj in ground truth to alignment method ...
            ret_value.append({'read_id1': query, 'read_id2': m})
            match_found += 1
    #print(match_found, len(overlapping_pair))
    os.remove(input_file)
    return ret_value


def process_read(header, sequence, ct):
    if header is None:
        return
    
    ### Run BLAST for the current read
    read_id = header.strip()
    read_seq = "".join(sequence).strip()

   # print(f'Creating temporary FASTA file for read: {ct}...')
    with open(f'read{ct}.fa', 'w') as temp_f:
        temp_f.write(f">{read_id}\n{read_seq}\n")

    output_file = f"blast_alignment_read{ct}.txt"
    
    blast_command = (
        f"blastn -db other_reads_db -query read{ct}.fa "
        f"-outfmt '6 qacc sacc sstart send' -out {output_file} -strand plus"
    )

    print(f"-> Running BLAST for read: {ct}")
    os.system(blast_command)

    os.remove(f'read{ct}.fa')
    ### convert to csv file
    rows = []

    with open(output_file, 'r') as f:
        for line in f:
            line = line.strip()


            fields = line.split()
            read1_id = fields[0]
            read2_id = fields[1]

            start = int(fields[2])
            end = int(fields[3])

            rows.append({
                'read1_id': read1_id,
                'read2_id': read2_id,
                'start': start,
                'end': end,
            })


    input_file = f'blast_alignment_read{ct}.csv'
    df = pd.DataFrame(rows)
    df.to_csv(input_file, sep='\t', index=False)
    os.remove(output_file)

    pairs = overlapping_pairs_for_read(input_file, 'overlaps_ground_truth.csv')
    all_overlapping_pairs.extend(pairs)
    df_pairs = pd.DataFrame(all_overlapping_pairs)
    df_pairs.to_csv('overlapping_pairs_output.csv', index=False)
    #print(f"Successfully saved all accuracy results ({len(df_pairs)} rows) to overlapping_pairs_output.csv")

### process each read from fasta file
def output_alignment_custom_files(fasta_file):
    read_ct = 0
    current_header = None
    current_sequence = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('>'):

                process_read(current_header, current_sequence, read_ct)
                read_ct += 1

                current_header = line[1:] 
                current_sequence = []
            else:
                current_sequence.append(line)

        process_read(current_header, current_sequence, read_ct)
        read_ct += 1
    print(f'Finished processing {read_ct} reads')
    

if __name__ == "__main__":
    fasta_file = 'readsMappingToChr1.fa'
    start_time = time.time()
    pairs = output_alignment_custom_files(fasta_file)
    end_time = time.time()
    total_time = end_time - start_time
    print(f"Total time taken: {total_time:.2f} seconds or {total_time/60:.2f} minutes")