from Bio import SeqIO
from lamp_coverage import primer_matches_sequence

part1 = "GCTGCGTTGTGTCTTGGGAGG"
part2 = "CTGTACGCATGGGGTAGC"

match_found = False
count = 0
for record in SeqIO.parse("sequences_DENGUE_1_NCBI_03_04_2026.fasta", "fasta"):
    t = str(record.seq).upper().replace('-', '')
    count += 1
    p1 = primer_matches_sequence(t, part1, max_errors=2, strict_3prime_len=0)
    p2 = primer_matches_sequence(t, part2, max_errors=2, strict_3prime_len=0)
    if p1 and p2:
        print("Found match on target %d" % count)
        match_found = True
        break

if not match_found:
    print("NO MATCH in entire DB!")
