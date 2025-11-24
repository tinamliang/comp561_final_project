1. Characterize the common types of sequencing error the pacbio machine makes using reads-to-entire-genome BLAST.
2. For each read Ri, run BLAST, where the database is the rest of the reads. 
3. Find all possible Rj that is the most similar to Ri.
4. Run a variant of S-W (ie. reduce weight of indel errors to substitution errors) algorithm to align the sequences and get a score. Get corrected overlap lengths.
5. Compare (Ri, Rj) results to the 'ground truth' obtained from BLAST.