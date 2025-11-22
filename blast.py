import pickle 

def read_file(fasta_file):
    genome = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            genome.append(line.strip())
    return ''.join(genome)

def read_db():
    # take read1 and run blast against database queried from pkl file.
    with open('index_w11.pkl', 'rb') as f:
        loaded_db = pickle.load(f)
    return loaded_db


def parse_chr1(fasta_file):
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            seq = line.strip()
            return seq
        
def blast(seq, db, genome):
    w = 11
    N = len(seq)
    for i in range(N-w+1):
        wmer = seq[i:i + w]
        if wmer in db:
            positions = db[wmer]
            #print(positions)
            for pos in positions:
                # sub_seq = seq[i:i + w]
                # sub_db_seq = genome[pos:pos + w]
                #print(f'q: {sub_seq}', f'D: {sub_db_seq}')
                l, r = findHSP(seq, genome, i, pos)
                #print(l, r)
     
def findHSP(seq, genome, i, pos):

    # get last indexes of the initial match
    new_i = i + 10
    new_pos = pos + 10

    left_extension, right_extension = True, True
    left_ind, right_ind = 0, 0
    max_left_score, max_right_score = 11, 11
    delta = 5
    left_score, right_score = 11, 11

    # final indices of the highest scoring pair
    f_l_ind, f_r_ind = i, new_i

    while left_extension or right_extension:

        if left_extension:
            print('extending left...')
            left_ind += 1

            # out of bounds
            if (i-left_ind) < 0 or (pos-left_ind) < 0:
                print('stopping left extension, index out of bound!')
                left_extension = False
                continue
            
            # match
            if seq[i-left_ind] == genome[pos-left_ind]:
                left_score += 1
            # mismatch
            else:
                left_score -= 2

            print(f'running sequence: {seq[i-left_ind:new_i+1], genome[pos-left_ind:new_pos+1]}')
            print(f'running max score: {max_left_score}, current score: {left_score}, f_l_ind: {f_l_ind}, left_ind: {left_ind}')
            
            # keep running max score
            if left_score >= max_left_score:
                max_left_score = left_score
                f_l_ind -= left_ind 

            # check delta threshold
            if left_score < max_left_score - delta:
                print(f'stopping left extension, delta reached! {f_l_ind}')
                left_extension = False

        if right_extension:
            print('extending right...')
            right_ind += 1

            # out of bounds
            if (new_i+right_ind) >= len(seq) or (new_pos+right_ind) >= len(genome):
                print('stopping right extension, index out of bound!')
                right_extension = False
                continue
            
            # match
            if seq[new_i+right_ind] == genome[new_pos+right_ind]:
                right_score += 1
            # mismatch
            elif seq[new_i+right_ind] != genome[new_pos+right_ind]:
                right_score -= 2
            
            print(f'running sequence: {seq[i:new_i+right_ind+1], genome[pos:new_pos+right_ind+1]}')
            print(f'running max score: {max_right_score}, current score: {right_score}, f_r_ind: {f_r_ind}, right_ind: {right_ind}')
            
            # keep running max score
            if right_score >= max_right_score:
                max_right_score = right_score
                f_r_ind += right_ind

            # check delta threshold
            if right_score < max_right_score - delta:
                print(f'stopping right extension, delta reached! {f_r_ind}')
                right_extension = False

    # l, r indices of the highest scoring pair.
    return f_l_ind, f_r_ind
  
if __name__ == "__main__":

    db_file = "GCA_000227135.2_ASM22713v2_genomic.fna"
    genome = read_file(db_file)

    fasta_file = "readsMappingToChr1.fa"
    db = read_db()
    seq = parse_chr1(fasta_file)
    blast(seq, db, genome)
