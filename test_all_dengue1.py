from Bio import SeqIO
from lamp_coverage import primer_matches_sequence

fip = "CATCCTGTCTGAAGCATTGGCTGGACAATTGACATGGAATGATC"
match_found = False
count = 0
for record in SeqIO.parse("sequences_DENGUE_1_NCBI_03_04_2026.fasta", "fasta"):
    t = str(record.seq).upper().replace('-', '')
    count += 1
    for allowed_err in [0, 1, 2]:
        for i in range(18, 26):
            part1 = fip[:i]
            part2 = fip[i:]
            p1 = primer_matches_sequence(t, part1, max_errors=allowed_err, strict_3prime_len=0)
            p2 = primer_matches_sequence(t, part2, max_errors=allowed_err, strict_3prime_len=0)
            if p1 and p2:
                print("Found match on target %d with err %d, split at i=%d" % (count, allowed_err, i))
                print("part1:", part1)
                print("part2:", part2)
                match_found = True
                break
        if match_found: break
    if match_found: break
