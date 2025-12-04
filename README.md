1. Obtain ground truth for read-to-read alignments using read-to-genome BLAST. 
2. Output csv file where each read Ri is identified with the Rj it overlaps with in the genome (overlaps_ground_truth.csv).
3. For each read Ri, run BLAST, where the database is the rest of the reads. 
4. Identify all high-scoring alignments. (results.csv)
5. Output a final csv identifying read Ri and the percentage of overlapping pairs the alignment method (Step 4) was able to correctly identify from the ones obtained in Step 2 (eval_naive_sol.csv).


### To try later ...
6. Run a variant of S-W (ie. reduce weight of indel errors to substitution errors) algorithm to align the sequences and get a score. Get corrected overlap lengths.
7. Compare (Ri, Rj) results to the 'ground truth' obtained from BLAST.
modifiedBLAST_pairing.py pairs highly aligned reads
