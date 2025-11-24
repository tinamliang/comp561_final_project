1. Characterize the common types of sequencing error the pacbio machine makes using reads-to-entire-genome BLAST.
2. For each read Ri, take its last k-mer and form a hash table with it. For each read Rj, determine
if there's a matching prefix k-mer. 
3. Confirm overlap and check that its length is above a minimum value.
4. Run a variant of S-W (ie. reduce weight of indel errors to substitution errors) algorithm to align the sequences and get a score. Get corrected overlap lengths.
5. Compare (Ri, Rj) results to the 'ground truth' obtained from BLAST.