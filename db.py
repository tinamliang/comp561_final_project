from collections import defaultdict
import pickle 

# from https://bioinformatics.stackexchange.com/questions/225/uppercase-vs-lowercase-letters-in-reference-genome: 
# 'N' are for hard-masked areas where the genome could not be assembled.
# lowercase letters are for soft-masked areas that are repetitive sequences.


# right now, upper to lower case regions are differentiated.
# do we want reads to align to highly likely intron regions? probably not.

def read_file(fasta_file):
    genome = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            genome.append(line.strip())
    return ''.join(genome)


def buildDB(genome, w=11):
    index_db = defaultdict(list)
    N = len(genome)
    for i in range(N-w+1):
        wmer = genome[i:i + w]
        if 'N' not in wmer:
            index_db[wmer].append(i)
    return index_db

def save_db(index_db, w=11):
    with open(f'index_w{w}.pkl', 'wb') as f:
        pickle.dump(index_db, f)
    print("Index saved.")


if __name__ == "__main__":
    fasta_file = "GCA_000227135.2_ASM22713v2_genomic.fna"
    genome = read_file(fasta_file)
    index_db = buildDB(genome, w=11)
    save_db(index_db)
