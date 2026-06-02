from Bio import SeqIO
lens = []
for r in SeqIO.parse("dengue_1_primers_aligned_26_05.fasta", "fasta"):
    l = len(str(r.seq).replace("-", ""))
    lens.append((r.id, l))
lens.sort(key=lambda x: x[1])
for i in range(20):
    print(lens[i])
