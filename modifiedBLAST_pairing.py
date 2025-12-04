import subprocess
import tempfile
# Insert sequencer reads in subjects
subject = "reads.fasta"

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
            header = line        
            seq_chunks = []
        else:
            seq_chunks.append(line)

    # save last record
    if header is not None:
        queries.append(header + "\n" + "".join(seq_chunks))

     
pairs = []
for i, q in enumerate(queries):
    with tempfile.NamedTemporaryFile("w+", suffix=".fasta", delete=False) as tmp:
        tmp.write(q + "\n")
        tmp.flush()
        blast_cmd = [
            "blastn",
            "-query", tmp.name,
            "-subject", subject,
            "-outfmt", "6",
            # "-use_sw_tback",
            "-reward", "1",
            "-penalty", "-3",
            "-gapopen", "2",
            "-gapextend", "1",
        ]

        result = subprocess.run(
            blast_cmd,
            capture_output=True,
            text=True,      # decode bytes -> str
        )                   # note: NO check=True here

        if not (len(result.stdout.strip().split('\n')) < 3):
            Ri = result.stdout.strip().split('\n')[2].split('\t')[0]
            Rj = result.stdout.strip().split('\n')[2].split('\t')[1]
            print((Ri, Rj))
            if Ri != Rj and reversed((Ri, Rj)) not in pairs:
                pairs.append((Ri, Rj))
print(len(pairs))

