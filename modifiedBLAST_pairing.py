import subprocess
import tempfile
import time

subject = r"reads.fasta"

# thresholds for calling an overlap from BLAST
MIN_LEN = 500      # min alignment length (tune this)
MIN_IDENT = 75.0   # min % identity (tune this for PacBio errors)

# -------------------------
# 1. Parse FASTA into queries
# -------------------------
queries = []
with open(subject, "r") as f:
    header = None
    seq_chunks = []

    for line in f:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            # save previous record
            if header is not None:
                queries.append(header + "\n" + "".join(seq_chunks))
            header = line            # keep '>'
            seq_chunks = []
        else:
            seq_chunks.append(line)

    # save last record
    if header is not None:
        queries.append(header + "\n" + "".join(seq_chunks))

print(f"Number of queries: {len(queries)}")

# -------------------------
# 2. Run BLAST per read and collect all overlapping pairs
# -------------------------
start_time = time.time()

pairs = set()  # store unique (id1, id2) with id1 < id2

for i, q in enumerate(queries):
    # write this read to a temp FASTA
    with tempfile.NamedTemporaryFile("w+", suffix=".fasta", delete=False) as tmp:
        tmp.write(q + "\n")
        tmp.flush()

        blast_cmd = [
            "blastn",
            "-query", tmp.name,
            "-subject", subject,
            "-outfmt", "6",        
            # scoring params to favour indels vs mismatches
            "-reward", "1",
            "-penalty", "-3",
            "-gapopen", "2",
            "-gapextend", "1",
            # optional speed tweaks:
            # "-dust", "no",
            # "-word_size", "11",
        ]

        result = subprocess.run(
            blast_cmd,
            capture_output=True,
            text=True,
        )

    # skip if BLAST crashed
    if result.returncode != 0:
        continue

    stdout = result.stdout.strip()
    if not stdout:
        continue

    for line in stdout.splitlines():
        fields = line.split('\t')
        if len(fields) < 12:
            continue

        qid = fields[0]
        sid = fields[1]
        pident = float(fields[2])
        length = int(fields[3])
        # other fields (qstart, qend, sstart, send, evalue, bitscore) are fields[6:12]

        if qid == sid:
            continue
        if length < MIN_LEN:
            continue
        if pident < MIN_IDENT:
            continue

        a, b = sorted((qid, sid))
        pairs.add((a, b))
        print(a, b, pident, length)

end_time = time.time()

print(f"Number of BLAST overlap pairs (unique): {len(pairs)}")
print(f"Execution time: {end_time - start_time:.2f} seconds")

# -------------------------
# 3. Write pairs to file (real read IDs)
# -------------------------
out_path = r"pairs_blast.csv"
with open(out_path, "w") as file_object:
    for a, b in sorted(pairs):
        file_object.write(f"{a},{b}\n")

