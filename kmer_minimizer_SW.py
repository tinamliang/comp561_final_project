import parasail
from Bio import SeqIO
from collections import Counter
import time

# Settings for parasail SW alignmnet 
MATCH = 1
MISMATCH = -3
GAP_OPEN = 2      # parasail uses positive gap_open, negative applied internally
GAP_EXTEND = 1

# Create a simple DNA scoring matrix equivalent to match/mismatch
matrix = parasail.matrix_create("ACGT", MATCH, MISMATCH)

def parasail_sw_score(q, s):
    """
    Local SW score using parasail's striped implementation with affine gaps.
    """
    if not q or not s:
        return 0

    result = parasail.sw_striped_16(
        q,
        s,
        GAP_OPEN,
        GAP_EXTEND,
        matrix
    )
    return result.score

# -----------------------------
# Parameters to tune
# -----------------------------
K = 15          # minimizer k-mer length
W = 30          # minimizer window size (in bases)
MAX_OCC = 200   # drop minimizers more frequent than this
MIN_SHARED = 5 # min shared minimizer hits to even consider SW
TOP_CANDS = 100  # max candidates per read (by shared minimizers)
SW_WINDOW = 2500   # suffix/prefix window size for SW
SCORE_THRESH = 120 # SW score cutoff for calling an overlap

FASTA_PATH = r"reads.fasta"
OUT_PATH   = r"pairs_kmer_SW.txt"


# -----------------------------
# Minimizers + index
# -----------------------------
def minimizers_in_read(seq, k, w):
    """Yield (minimizer, position) for each window."""
    assert w >= k
    n = len(seq)
    if n < k or n < w:
        return

    kmers = [seq[i:i + k] for i in range(0, n - k + 1)]

    # number of windows = n - w + 1
    for start in range(0, len(kmers) - (w - k) + 1):
        window = kmers[start:start + (w - k + 1)]
        # choose lexicographically smallest k-mer
        min_kmer = min(window)
        rel_pos = window.index(min_kmer)
        pos_in_read = start + rel_pos
        yield min_kmer, pos_in_read


def build_minimizer_index(fasta_path, k, w, max_occurrences=None):
    """Build minimizer index: minimizer -> list of (read_idx, position)."""
    kmer_index = {}
    reads = []

    for i, rec in enumerate(SeqIO.parse(fasta_path, "fasta")):
        seq = str(rec.seq).upper()
        reads.append((rec.id, seq))

        for kmer, pos in minimizers_in_read(seq, k, w):
            if "N" in kmer:
                continue
            kmer_index.setdefault(kmer, []).append((i, pos))

    # Optional: drop super-common minimizers (repeats)
    if max_occurrences is not None:
        kmer_index = {
            kmer: hits
            for kmer, hits in kmer_index.items()
            if len(hits) <= max_occurrences
        }

    return reads, kmer_index


def find_overlap_candidates(read_id, reads, kmer_index, k, w):
    """
    Use minimizers on the query read to find candidate overlapping reads.
    Return (candidates_list, counts_by_read_idx).
    """
    _, seq = reads[read_id]
    counts = Counter()

    for kmer, pos in minimizers_in_read(seq, k, w):
        if "N" in kmer:
            continue
        if kmer not in kmer_index:
            continue
        for other_id, other_pos in kmer_index[kmer]:
            if other_id == read_id:
                continue
            counts[other_id] += 1

    # we don't threshold here yet; caller decides MIN_SHARED
    candidates = list(counts.keys())
    return candidates, counts


# -----------------------------
# Overlap scoring (SW on suffix/prefix)
# -----------------------------
dWINDOW = 2500   # or 2000, tune this

def get_suffix_prefix_window(qseq, sseq, window=SW_WINDOW):
    q_sub = qseq[-window:] if len(qseq) > window else qseq
    s_sub = sseq[:window]  if len(sseq) > window else sseq
    return q_sub, s_sub


# -----------------------------
# Main
# -----------------------------
start_time = time.time()

# Build index 
reads, kmer_index = build_minimizer_index(
    FASTA_PATH, k=K, w=W, max_occurrences=MAX_OCC
)


# For each read, query the existing index and run SW on good candidates
overlap_pairs = set()     # store (read_id_str1, read_id_str2)

for read_id in range(len(reads)):
    qname, qseq = reads[read_id]

    # Find candidates via shared minimizers
    candidates, counts = find_overlap_candidates(
        read_id, reads, kmer_index, k=K, w=W
    )

    # Keep only strong candidates (enough shared minimizers), limited to TOP_CANDS
    strong = [
        (rid, c) for rid, c in counts.most_common(TOP_CANDS)
        if c >= MIN_SHARED
    ]
    if not strong:
        continue

    for cand_id, c in strong:
        sname, sseq = reads[cand_id]

        # restrict SW to suffix(Ri) vs prefix(Rj)
        q_sub, s_sub = get_suffix_prefix_window(qseq, sseq, window=SW_WINDOW)
        score = parasail_sw_score(q_sub, s_sub)

        if score >= SCORE_THRESH:
            a, b = sorted((qname, sname))
            if (a, b) not in overlap_pairs:
                overlap_pairs.add((a, b))
            print(a, b, score)


end_time = time.time()
print(f"Execution time: {end_time - start_time:.2f} seconds")
print(f"Number of overlap pairs (unique): {len(overlap_pairs)}")

# Write pairs to file using real read IDs
with open(OUT_PATH, "w") as f:
    for a, b in sorted(overlap_pairs):
        f.write(f"{a},{b}\n")
