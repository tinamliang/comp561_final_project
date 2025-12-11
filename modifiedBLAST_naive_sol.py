import os
from collections import defaultdict
import pandas as pd

all_accuracy_results = []
def calculate_accuracy_for_read(input_file, output_file):
    query = None
    matches = defaultdict(list)
    match_found = 0
    overlapping_pair = []

    with open(input_file, 'r') as f:
        next(f)
        next(f)
        for line in f:
            line = line.strip()

            # Split by any whitespace
            fields = line.split()
            query = fields[0]
            if fields[1] not in matches[query]:
                matches[query].append(fields[1])
            
   # print(f'query: {query}, matches: {matches}')

    with open(output_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Split by any whitespace
            fields = line.split()
            q = fields[0]

            if query == q:
                overlapping_pair.append(fields[1])
    

    for m in matches[query]:
        if m in overlapping_pair:
            match_found += 1
        #   print(f'Match found: {m}')
    
    if overlapping_pair:
        accuracy = match_found / len(overlapping_pair)
        print(f'Accuracy= {match_found}/{len(overlapping_pair)}: {(accuracy * 100):.2f}%')
        os.remove(input_file)
        return {
            'query': query,
            'accuracy': f'{(accuracy * 100):.2f}%'
        }
    else:
        os.remove(input_file)
        return {
            'query': query,
            'accuracy': '0%'
        }


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
        f"-outfmt '6 qacc sacc sstart send' -out {output_file}"
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
   # print(f"Successfully parsed {len(df)} lines and saved to {input_file}")
    os.remove(output_file)
  #  print(f'removing temporary file: {output_file}')

    read_accuracy_data = calculate_accuracy_for_read(input_file, 'overlaps_ground_truth.csv')
    all_accuracy_results.append(read_accuracy_data)

def accuracy_table():
    if all_accuracy_results:
        #print(all_accuracy_results)
        df_accuracy = pd.DataFrame(all_accuracy_results)
        df_accuracy.to_csv('results_with_both_strands.csv', index=False)
        print(f"Successfully saved all accuracy results ({len(df_accuracy)} rows) to results_with_both_strands.csv")

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
    output_alignment_custom_files(fasta_file)
    accuracy_table()