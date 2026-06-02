from Bio import SeqIO
from lamp_coverage import primer_matches_sequence

fip = "TTGGGCCCCCATTGTTGCTGTTTTAGTGGACTAGCGGTTAGAGG"
targets = []
for record in SeqIO.parse("sequences_DENGUE_2_NCBI_03_04_2026.fasta", "fasta"):
    t = str(record.seq).upper().replace('-', '')
    if len(t) > 100:
        targets.append(t)
        if len(targets) >= 50:
            break

found = False
for count, t in enumerate(targets):
    if found: break
    for allowed_err in [0, 1, 2]:
        if found: break
        for linker_len in [0, 4, 1, 2, 3, 5, 6, 7, 8, 9]:
            if found: break
            for i in range(14, len(fip)-13):
                if i + linker_len > len(fip) - 14:
                    continue
                if linker_len > 0:
                    linker_seq = fip[i:i+linker_len]
                    if len(set(linker_seq)) != 1:
                        continue
                p1 = fip[:i]
                p2 = fip[i+linker_len:]
                pos1 = primer_matches_sequence(t, p1, max_errors=allowed_err, strict_3prime_len=0)
                pos2 = primer_matches_sequence(t, p2, max_errors=allowed_err, strict_3prime_len=0)
                if pos1 and pos2:
                    print(f"Match on target {count} | err={allowed_err} | linker={linker_len} | i={i}")
                    print(f"p1: {p1}")
                    print(f"p2: {p2}")
                    found = True
                    break
