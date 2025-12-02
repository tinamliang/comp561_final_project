### some script code was produced using Claude: with prompt "can we parse just the gaps metrics and identities % from txt file"
"""
makeblastdb -in GCA_000227135.2_ASM22713v2_genomic.fna -dbtype nucl -out leishmania_db
blastn -query readsMappingToChr1.fa -db leishmania_db -strand plus -out alignments.txt -outfmt 0
   alignments.txt: full read-to-genome BLAST alignment output

blastn -query readsMappingToChr1.fa -db leishmania_db -strand plus -out alignments_custom_tab.txt -outfmt "6 qseqid sstart send"
  alignments_custom_tab.txt: custom tabular output with read start/end positions on genome

"""
import re

def parse_blast_metrics(blast_text):
    """Parse Identities and Gaps from BLAST output."""
    results = []
    
    # Split by Query= to get individual queries
    queries = re.split(r'Query=\s*', blast_text)
    
    for query in queries[1:]:  # Skip first empty split
        # Get query name (first line)
        query_name = query.split('\n')[0].strip()
        
        # Find all Identities/Gaps lines
        identity_pattern = r'Identities\s*=\s*(\d+)/(\d+)\s*\((\d+)%\)'
        gaps_pattern = r'Gaps\s*=\s*(\d+)/(\d+)\s*\((\d+)%\)'
        
        identities = re.findall(identity_pattern, query)
        gaps = re.findall(gaps_pattern, query)

        for i, (ident, gap) in enumerate(zip(identities, gaps)):
            results.append({
                'query': query_name[:50] + '...' if len(query_name) > 50 else query_name,
                'identity_matches': f"{ident[0]}/{ident[1]}",
                'identity_pct': f"{ident[2]}%",
                'gaps_count': f"{gap[0]}/{gap[1]}",
                'gaps_pct': f"{gap[2]}%"
            })
    
    return results

if __name__ == "__main__":

    with open('alignments.txt', 'r') as f:
        blast_output = f.read()

    metrics = parse_blast_metrics(blast_output)
    print(len(metrics))
    
    # with open('metrics.txt', 'w') as f:
    #     for m in metrics:
    #         print(f"Query: {m['query']}", file=f)
    #         print(f"  Matches (%): {m['identity_pct']}", file=f)
    #         print(f"  Gaps (%): {m['gaps_pct']}", file=f)
    #         subst_pct = f"{100 - int(m['identity_pct'][:-1]) - int(m['gaps_pct'][:-1])}%"
    #         print(f"  Substitutions (%): {subst_pct}", file=f)
    #         print()