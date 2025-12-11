1. Obtain ground truth for read-to-read alignments using read-to-genome BLAST. 
2. Output csv file where each read Ri is identified with the Rj it overlaps with in the genome (overlaps_ground_truth.csv).
3. For each read Ri, run BLAST, where the database is the rest of the reads. 
4. Identify all high-scoring alignments. (results.csv)
5. Output a final csv identifying read Ri and the percentage of overlapping pairs the alignment method (Step 4) was able to correctly identify from the ones obtained in Step 2 (eval_naive_sol.csv).
6. Run a modified BLAST (ie. reduce weight of indel errors to substitution errors) algorithm to align the sequences and get a score. Choose pairs based on minimum overlap length and % identity.
7. Run a K-mer/minimizer with modified Smith-Waterman (same parameters as modified BLAST). Choose pairs based on candidates filtered out by minimizer and kmer indexing and above SW score threshold.  
8. Compare (Ri, Rj) from the different approaches to the 'ground truth' obtained from BLAST.

