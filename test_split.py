import sys
from Bio import SeqIO
from collections import defaultdict
from lamp_coverage import primer_matches_sequence, seq_to_regex, IUPAC_DICT, auto_split_fip_bip

targets = {}
for record in SeqIO.parse("sequences_DENGUE_1_NCBI_03_04_2026.fasta", "fasta"):
    clean_seq = str(record.seq).upper().replace('-', '')
    if len(clean_seq) > 100:
        targets[record.description] = clean_seq
        if len(targets) >= 5:
            break

primer_sets = {
    '1': {
        'FIP': 'TCCTCTGGCCAATGATATCTGCAGTRAACATGACATCYAGAA',
        'BIP': 'TGAACAYAAGTCAACATGGCACARTYTCACCACDCCATT'
    }
}
txt = {
    'split_success': "Auto-split réussi pour Set {} ({}) -> {} et {} (linker de {} nt ignoré).",
    'split_fail': "Attention : Impossible de déterminer la coupe automatique pour Set {} ({})."
}

auto_split_fip_bip(primer_sets, targets, txt)
print(primer_sets)
