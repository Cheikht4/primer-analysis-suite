#!/usr/bin/env python3
import argparse
import sys
import os
import itertools
import regex
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

# Import optionnel de tqdm pour les barres de progression / Optional tqdm import for progress bars
try:
    from tqdm import tqdm
except ImportError:
    # Fallback si tqdm n'est pas installé / Fallback if tqdm is not installed
    def tqdm(iterable=None, **kwargs):
        desc = kwargs.get('desc', '')
        if desc and iterable is None:
            class _DummyBar:
                def update(self, n=1): pass
                def set_postfix_str(self, s): pass
                def close(self): pass
            return _DummyBar()
        if desc:
            print(f"{desc}...")
        return iterable if iterable is not None else []

# Dictionnaire IUPAC vers expression régulière
IUPAC_DICT = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

def seq_to_regex(seq):
    """Convertit une séquence avec codes IUPAC en expression régulière."""
    pattern = ""
    for char in seq.upper():
        if char.isalpha():
            pattern += IUPAC_DICT.get(char, char)
    return pattern

def primer_matches_sequence(target_seq, primer_seq, max_errors, strict_3prime_len=3, strict_3prime_tolerate=0):
    """
    Vérifie si l'amorce match la séquence cible (brin sens ou anti-sens).
    Retourne un tuple (start, end) de la position du match sur le brin sens,
    ou None s'il n'y a pas de match.
    """
    # Définition des positions tolérées / Definition of tolerated positions
    # Position 1 = base -1, Position 2 = base -2, Position 3 = base -3 (compté depuis l'extrémité 3' / counted from 3' end)
    tolerate_positions = set()
    if strict_3prime_tolerate == 1:
        tolerate_positions = {2}
    elif strict_3prime_tolerate == 2:
        tolerate_positions = {1, 2}

    def build_regex(seq, max_e, strict_len, is_rc=False):
        if strict_len > 0 and len(seq) > strict_len:
            # Choix du type d'erreur pour la partie 5' : substitutions uniquement (s) en mode ARMS pour éviter les décalages d'indels
            # Choice of error type for the 5' part: substitutions only (s) in ARMS mode to avoid indel shifts
            type_e = 's' if strict_3prime_tolerate > 0 else 'e'
            
            if not is_rc:
                # Brin sens / Sense strand
                seq_5 = seq[:-strict_len]
                pattern_5 = seq_to_regex(seq_5)
                
                # Construction de la zone 3' base par base / Building the 3' region base by base
                pattern_3 = ""
                # De la base la plus éloignée de 3' vers la base terminale (5' vers 3') / From the base furthest from 3' to the terminal base (5' to 3')
                for k in range(strict_len, 0, -1):
                    base = seq[-k]
                    base_regex = seq_to_regex(base)
                    if k in tolerate_positions:
                        pattern_3 += f"(?:{base_regex}){{s<=1}}"
                    else:
                        pattern_3 += base_regex
                
                return f"(?e)(?:{pattern_5}){{{type_e}<={max_e}}}{pattern_3}"
            else:
                # Brin anti-sens / Antisense strand
                primer_rc = seq # Ici seq est déjà primer_rc / Here seq is already primer_rc
                seq_5_rc = primer_rc[strict_len:]
                pattern_5_rc = seq_to_regex(seq_5_rc)
                
                # Le 3' de l'amorce originale correspond au début (5') de primer_rc / The 3' of the original primer corresponds to the beginning (5') of primer_rc
                pattern_3_rc = ""
                for k in range(1, strict_len + 1):
                    base = primer_rc[k-1]
                    base_regex = seq_to_regex(base)
                    if k in tolerate_positions:
                        pattern_3_rc += f"(?:{base_regex}){{s<=1}}"
                    else:
                        pattern_3_rc += base_regex
                
                return f"(?e){pattern_3_rc}(?:{pattern_5_rc}){{{type_e}<={max_e}}}"
        else:
            pattern = seq_to_regex(seq)
            return f"(?e)({pattern}){{e<={max_e}}}"

    # Brin sens
    regex_pattern = build_regex(primer_seq, max_errors, strict_3prime_len, is_rc=False)
    match_sense = regex.search(regex_pattern, target_seq, regex.BESTMATCH)
    if match_sense:
        return match_sense.span()
        
    # Brin anti-sens
    primer_rc = str(Seq(primer_seq).reverse_complement())
    regex_rc_pattern = build_regex(primer_rc, max_errors, strict_3prime_len, is_rc=True)
    match_antisense = regex.search(regex_rc_pattern, target_seq, regex.BESTMATCH)
    
    if match_antisense:
        return match_antisense.span()
        
    return None

def load_primers(filepath, is_pcr=False):
    """
    Charge les amorces et les regroupe par set.
    Attend un format de nom du type >SetID_PrimerName (ex: >2_F3).
    """
    primer_sets = defaultdict(dict)
    
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline()
            if not first_line: return {}
            
            f.seek(0)
            first_char = f.read(1)
            f.seek(0)
            
            records = []
            if first_char == '>':
                records = list(SeqIO.parse(f, "fasta"))
            else:
                # Fichier texte basique
                primer_count = 1
                for line in f:
                    line = line.strip()
                    if not line: continue
                    parts = line.split()
                    if len(parts) >= 2:
                        name = parts[0]
                        seq = "".join(parts[1:])
                    else:
                        name = f"DefaultSet_Primer{primer_count}"
                        seq = parts[0]
                        primer_count += 1
                    
                    seq = "".join(c for c in seq if c.isalpha())
                    records.append({'id': name, 'seq': seq})
                    
            for record in records:
                if hasattr(record, 'id'):
                    name = record.id
                    seq = str(record.seq).upper()
                else:
                    name = record['id']
                    seq = record['seq'].upper()
                    
                if '_' in name:
                    # Découpage sur le DERNIER underscore pour associer les amorces entre elles
                    # Split on the LAST underscore to group primers by set
                    # Ex: "SetA_Dengue2_F3" → set_id="SetA_Dengue2", primer_id="F3"
                    parts = name.rsplit('_', 1)
                    set_id = parts[0]
                    primer_id = parts[1]
                else:
                    set_id = "Default"
                    primer_id = name
                    
                # Gestion des alias
                if is_pcr:
                    alias_map = {
                        'FWD': 'F', 'FORWARD': 'F', 'F3': 'F', 'FP': 'F',
                        'REV': 'R', 'REVERSE': 'R', 'B3': 'R', 'RP': 'R',
                        'PROBE': 'P', 'SONDE': 'P'
                    }
                else:
                    alias_map = {
                        'BLP': 'BLOOP',
                        'FLP': 'FLOOP',
                        'LOOPF': 'FLOOP',
                        'LOOPB': 'BLOOP',
                        'LF': 'FLOOP',
                        'LB': 'BLOOP',
                        'F1C': 'F1',
                        'B1C': 'B1'
                    }
                primer_id = alias_map.get(primer_id.upper(), primer_id)
                    
                primer_sets[set_id][primer_id] = seq
                
    except Exception as e:
        print(f"Error reading primers file / Erreur lecture fichier amorces : {e}")
        sys.exit(1)
        
    return primer_sets

def auto_split_fip_bip(primer_sets, targets, txt):
    """
    Tente de séparer automatiquement les amorces nommées FIP ou BIP
    en trouvant l'alignement de leurs sous-parties sur les premières séquences cibles.
    """
    # Prendre jusqu'à 200 cibles pour tester (au cas où les premières auraient des mutations divergentes)
    test_targets = list(targets.values())[:200]
    if not test_targets: return
    
    for set_id, primers in list(primer_sets.items()):
        for combo_name in ['FIP', 'BIP']:
            if combo_name in primers:
                seq_combo = primers[combo_name]
                found_split = False
                
                # On cherche d'abord avec 0 erreur, puis 1, puis 2
                for allowed_err in [0, 1, 2]:
                    if found_split: break
                    # Ensuite on teste l'absence de linker (0), puis les linkers biologiquement valides (4 à 9 bases)
                    # Next we test the absence of linker (0), then biologically valid linkers (4 to 9 bases)
                    for linker_len in [0, 4, 5, 6, 7, 8, 9]:
                        if found_split: break
                        for t_seq in test_targets:
                            if found_split: break
                            
                            valid_splits = []
                            L = len(seq_combo)
                            for i in range(14, L - 13):
                                if i + linker_len > L - 14:
                                    continue
                                
                                if linker_len > 0:
                                    linker_seq = seq_combo[i:i+linker_len]
                                    # Contrainte biologique: un linker doit être une répétition de base (ex: TTTT, AAAA)
                                    if len(set(linker_seq)) != 1:
                                        continue
                                
                                part1 = seq_combo[:i]
                                part2 = seq_combo[i+linker_len:]
                                
                                pos1 = primer_matches_sequence(t_seq, part1, max_errors=allowed_err, strict_3prime_len=0)
                                pos2 = primer_matches_sequence(t_seq, part2, max_errors=allowed_err, strict_3prime_len=0)
                                
                                if pos1 and pos2:
                                    dist = abs(pos2[0] - pos1[0])
                                    if dist < 250:
                                        if combo_name == 'FIP' and pos2[0] < pos1[0]:
                                            valid_splits.append((part1, part2))
                                        elif combo_name == 'BIP' and pos1[0] < pos2[0]:
                                            valid_splits.append((part1, part2))
                        
                            if valid_splits:
                                # En cas d'ambiguïté (plusieurs coupes valides avec 0 erreur), 
                                # on choisit la coupe qui équilibre le mieux la taille des deux amorces.
                                best_split = min(valid_splits, key=lambda x: abs(len(x[0]) - len(x[1])))
                                part1, part2 = best_split
                                
                                if combo_name == 'FIP':
                                    primers['F1'] = str(Seq(part1).reverse_complement())
                                    primers['F2'] = part2
                                else:
                                    primers['B1'] = str(Seq(part1).reverse_complement())
                                    primers['B2'] = part2
                                    
                                found_split = True
                                print(txt['split_success'].format(set_id, combo_name, combo_name[0]+'1', combo_name[0]+'2', linker_len))
                
                if not found_split:
                    import re
                    # Fallback heuristique : si la souche compatible n'est pas dans les 200 premières,
                    # on tente une coupe "à l'aveugle" en cherchant au moins 4 "T" au milieu (méthode manuelle classique)
                    match = re.search(r'T{4,}', seq_combo[14:-14])
                    if match:
                        linker_start = match.start() + 14
                        part1 = seq_combo[:linker_start]
                        part2 = seq_combo[linker_start+4:]  # On ignore exactement 4 "T"
                        
                        if combo_name == 'FIP':
                            primers['F1'] = str(Seq(part1).reverse_complement())
                            primers['F2'] = part2
                        else:
                            primers['B1'] = str(Seq(part1).reverse_complement())
                            primers['B2'] = part2
                            
                        found_split = True
                        print(txt['split_success'].format(set_id, combo_name, combo_name[0]+'1', combo_name[0]+'2', '4 (fallback heuristique)'))
                
                if found_split:
                    del primers[combo_name]
                else:
                    print(txt['split_fail'].format(set_id, combo_name))

def main():
    parser = argparse.ArgumentParser(description="Évalue la couverture et l'ordre des amorces LAMP / Evaluate LAMP primer coverage and structural order.")
    parser.add_argument("-t", "--target", required=True, help="Fichier FASTA cible / Target FASTA file.")
    parser.add_argument("-p", "--primers", required=True, help="Fichier FASTA ou TXT des amorces / Primers FASTA or TXT file.")
    parser.add_argument("-o", "--output", required=True, help="Fichier de rapport en sortie / Output report file.")
    parser.add_argument("-e", "--errors", type=int, default=0, help="Nombre max d'erreurs hors zone 3' / Max errors outside 3' region. Def: 0")
    parser.add_argument("-s", "--strict-3prime", type=int, default=3, dest="strict_3prime", help="Taille zone 3' stricte / Strict 3' region size. Def: 3")
    parser.add_argument("--strict-3prime-tolerate", type=int, choices=[0, 1, 2], default=0, help="Niveau de tolérance en zone 3' (0: tout strict, 1: pos 2 tolérée, 2: pos 1 et 2 tolérées). / Tolerance level in the 3' region (0: all strict, 1: pos 2 tolerated, 2: pos 1 and 2 tolerated).")
    parser.add_argument("--strict-intersection", action="store_true", help="Exige que toutes les amorces du fichier matchent la cible (comportement strict historique). / Requires all primers in the file to match the target (historical strict behavior).")
    
    # Options de sortie
    parser.add_argument("--summary-only", action="store_true", help="N'affiche que les statistiques / Output only summary statistics.")
    parser.add_argument("--combine", action="store_true", help="Couverture combinatoire 2 à 2 / Calculate 2-by-2 multiplexing coverage.")
    parser.add_argument("--export-seqs", action="store_true", help="Exporte les séquences validées / Export validated sequences per set.")
    parser.add_argument("--pcr", action="store_true", help="Mode PCR : Gère les amorces Forward, Reverse et Sonde / PCR mode: handles Fwd, Rev and Probe primers.")
    
    # Langue
    parser.add_argument("--lng", type=str, default="en", choices=["en", "fr"], help="Langue du rapport / Report language (en, fr). Def: en")
    
    args = parser.parse_args()
    
    # Dictionnaire de traduction
    lang = args.lng.lower()
    T = {
        'fr': {
            'target_load': "Chargement des séquences cibles...",
            'target_err': "Erreur : Aucune séquence cible trouvée.",
            'target_loaded': "séquences cibles chargées.",
            'primer_load': "Chargement des amorces...",
            'primer_err': "Erreur : Aucune amorce trouvée.",
            'primer_loaded': "set(s) d'amorces identifié(s).",
            'analysis': "Analyse en cours avec tolérance = {} erreur(s) (et {} nt stricts en 3')...",
            'progression': "Progression : {}%",
            'report_gen': "Génération du rapport...",
            'report_title': "Rapport de Couverture des Amorces LAMP",
            'target_file': "Fichier cible",
            'primer_file': "Fichier amorces",
            'tolerance': "Tolérance d'erreurs par amorce",
            'strict_3': "avec {} nt stricts en 3'",
            'active_options': "Options actives",
            'set_title': "--- SET D'AMORCES : {} ---",
            'indiv_match': "Matchs individuels par amorce :",
            'global_raw': "Match Global du Set (Intersection Brute, toutes amorces présentes)",
            'global_valid': "Match Global Valide (Intersection + Ordre Correct structurel LAMP)",
            'amplified_seqs': "Séquences amplifiées théoriquement par le Set {} :",
            'table_header': "Séquence_ID\tTaille_Amplicon\tStatut_Ordre\tOrdre_Observe",
            'order_correct': "Ordre Correct",
            'order_incorrect': "Ordre Incorrect",
            'no_seq': "(Aucune séquence ne remplit les critères d'amplification)",
            'export_success': "-> Fichier d'export généré avec les séquences valides : {}",
            'export_err': "-> Erreur lors de l'export des séquences : {}",
            'combine_title': "Calcul Combinatoire (Multiplexage 2 à 2)",
            'combine_subtitle': "Couverture théorique si les sets sont combinés dans la même réaction (Union des séquences valides) :",
            'done': "Terminé. Rapport principal sauvegardé dans : {}",
            'write_err': "Erreur lors de l'écriture du rapport : {}",
            'seqs_word': "séquences",
            'split_success': "Auto-split réussi pour Set {} ({}) -> {} et {} (linker de {} nt ignoré).",
            'split_fail': "Attention : Impossible de déterminer la coupe automatique pour Set {} ({}).",
        },
        'en': {
            'target_load': "Loading target sequences...",
            'target_err': "Error: No target sequences found.",
            'target_loaded': "target sequences loaded.",
            'primer_load': "Loading primers...",
            'primer_err': "Error: No primers found.",
            'primer_loaded': "primer set(s) identified.",
            'analysis': "Analysis in progress with tolerance = {} error(s) (and {} strict 3' nt)...",
            'progression': "Progress : {}%",
            'report_gen': "Generating report...",
            'report_title': "LAMP Primer Coverage Report",
            'target_file': "Target file",
            'primer_file': "Primer file",
            'tolerance': "Error tolerance per primer",
            'strict_3': "with {} strict 3' nt",
            'active_options': "Active options",
            'set_title': "--- PRIMER SET : {} ---",
            'indiv_match': "Individual matches per primer:",
            'global_raw': "Set Global Match (Raw Intersection, all primers present)",
            'global_valid': "Set Valid Global Match (Intersection + Structurally Correct LAMP Order)",
            'amplified_seqs': "Theoretically amplified sequences by Set {} :",
            'table_header': "Sequence_ID\tAmplicon_Size\tOrder_Status\tObserved_Order",
            'order_correct': "Correct Order",
            'order_incorrect': "Incorrect Order",
            'no_seq': "(No sequences meet the amplification criteria)",
            'export_success': "-> Export file generated with valid sequences: {}",
            'export_err': "-> Error exporting sequences: {}",
            'combine_title': "Combinatorial Calculation (2-by-2 Multiplexing)",
            'combine_subtitle': "Theoretical coverage if sets are combined in the same reaction (Union of valid sequences):",
            'done': "Done. Main report saved in: {}",
            'write_err': "Error writing report: {}",
            'seqs_word': "sequences",
            'split_success': "Auto-split successful for Set {} ({}) -> {} and {} ({} nt linker ignored).",
            'split_fail': "Warning: Could not automatically determine split for Set {} ({}).",
        }
    }
    
    txt = T[lang]
    if args.pcr:
        if lang == 'fr':
            txt['report_title'] = "Rapport de Couverture des Amorces PCR"
            txt['global_valid'] = "Match Global Valide (Intersection + Ordre PCR Correct)"
        else:
            txt['report_title'] = "PCR Primer Coverage Report"
            txt['global_valid'] = "Set Valid Global Match (Intersection + Structurally Correct PCR Order)"
    
    print(txt['target_load'])
    targets = {}
    try:
        # Chargement avec barre de progression / Loading with progress bar
        all_records = list(tqdm(
            SeqIO.parse(args.target, "fasta"),
            desc="📂 Chargement / Loading" if args.lng == 'fr' else "📂 Loading sequences",
            unit=" seq",
            colour="cyan"
        ))
        for record in all_records:
            clean_seq = str(record.seq).upper().replace('-', '')
            if len(clean_seq) > 100:
                targets[record.description] = clean_seq
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    total_targets = len(targets)
    if total_targets == 0:
        print(txt['target_err'])
        sys.exit(1)
    print(f"{total_targets} {txt['target_loaded']}")
    
    print(txt['primer_load'])
    primer_sets = load_primers(args.primers, is_pcr=args.pcr)
    if not primer_sets:
        print(txt['primer_err'])
        sys.exit(1)
    
    suspected_pcr = False
    if not args.pcr:
        # Tentative de séparation auto de FIP et BIP
        auto_split_fip_bip(primer_sets, targets, txt)
        
        # Vérification si un essai PCR est suspecté en mode LAMP / Checking if a PCR assay is suspected in LAMP mode
        for s_id, primers in primer_sets.items():
            set_keys = {k.upper() for k in primers.keys()}
            if 'FP' in set_keys or 'RP' in set_keys or 'PROBE' in set_keys or (set_keys.issubset({'F', 'R', 'P'}) and not set_keys.intersection({'F3', 'B3', 'F2', 'F1', 'B1', 'B2', 'FIP', 'BIP'})):
                suspected_pcr = True
                break
                
        if suspected_pcr:
            print("\n" + "!" * 70)
            if args.lng == 'fr':
                print("⚠️  ATTENTION : Les amorces ressemblent à un essai PCR (FP, RP, Probe).")
                print("   Le script s'exécute actuellement en mode LAMP.")
                print("   Veuillez relancer avec l'option '--pcr' pour une analyse correcte.")
            else:
                print("⚠️  WARNING: The primers look like a PCR assay (FP, RP, Probe).")
                print("   The script is currently running in LAMP mode.")
                print("   Please rerun with the '--pcr' option for a correct analysis.")
            print("!" * 70 + "\n")
    
    print(f"{len(primer_sets)} {txt['primer_loaded']}")
    print(txt['analysis'].format(args.errors, args.strict_3prime))
    
    primer_matches = defaultdict(lambda: defaultdict(set))
    primer_positions = defaultdict(lambda: defaultdict(dict))

    # Boucle principale d'analyse avec barre de progression tqdm
    # Main analysis loop with tqdm progress bar
    total_steps = len(primer_sets) * len(targets)
    bar_label = "🔬 Analyse" if args.lng == 'fr' else "🔬 Analysing"

    with tqdm(total=total_steps, desc=bar_label, unit=" seq", colour="green") as pbar:
        for set_id, primers in primer_sets.items():
            for seq_id, seq in targets.items():
                for primer_id, primer_seq in primers.items():
                    pos = primer_matches_sequence(seq, primer_seq, args.errors, args.strict_3prime, args.strict_3prime_tolerate)
                    if pos:
                        primer_matches[set_id][primer_id].add(seq_id)
                        primer_positions[set_id][seq_id][primer_id] = pos

                pbar.update(1)
                # Affiche le set en cours dans la barre / Show current set in bar
                pbar.set_postfix_str(f"Set {set_id}")

    print(txt['report_gen'])
    
    if args.pcr:
        MASTER_ORDER = ['F', 'P', 'R']
    else:
        # Ordre général LAMP incluant à la fois les LOOP et les STEM (s'ils sont là)
        MASTER_ORDER = ['F3', 'F2', 'FLOOP', 'F1', 'STEMF', 'STEMB', 'B1', 'BLOOP', 'B2', 'B3']
    MASTER_ORDER_RC = list(reversed(MASTER_ORDER))
    
    # Définition des amorces essentielles (utilisées pour l'intersection relaxée)
    # Definition of essential primers (used for relaxed intersection)
    ESSENTIAL_LAMP = {'F3', 'B3', 'F2', 'F1', 'B1', 'B2', 'FIP', 'BIP'}
    ESSENTIAL_PCR = {'F', 'R'}
    
    # Stockage des sets de séquences valides pour le calcul combinatoire
    valid_sequences_per_set = {}
    
    try:
        with open(args.output, 'w') as out:
            if not args.pcr and suspected_pcr:
                if args.lng == 'fr':
                    out.write("⚠️  ATTENTION : Détecté comme essai PCR mais exécuté en mode LAMP. Veuillez utiliser --pcr.\n")
                else:
                    out.write("⚠️  WARNING: Detected as PCR assay but executed in LAMP mode. Please use --pcr.\n")
                out.write("=" * 80 + "\n\n")
                
            out.write(f"{txt['report_title']}\n")
            out.write("=" * max(len(txt['report_title']), 38) + "\n")
            out.write(f"{txt['target_file']} : {args.target} ({total_targets} {txt['seqs_word']})\n")
            out.write(f"{txt['primer_file']} : {args.primers}\n")
            out.write(f"{txt['tolerance']} : {args.errors} ({txt['strict_3'].format(args.strict_3prime)})\n")
            
            # Afficher les options actives
            active_options = []
            if args.summary_only: active_options.append("--summary-only")
            if args.combine: active_options.append("--combine")
            if args.export_seqs: active_options.append("--export-seqs")
            if args.strict_intersection: active_options.append("--strict-intersection")
            if active_options:
                out.write(f"{txt['active_options']} : {', '.join(active_options)}\n")
            out.write("\n")
            
            for set_id, primers in primer_sets.items():
                out.write(txt['set_title'].format(set_id) + "\n")
                
                set_matches_list = []
                out.write(f"{txt['indiv_match']}\n")
                for primer_id, primer_seq in primers.items():
                    matches = primer_matches[set_id][primer_id]
                    set_matches_list.append(matches)
                    match_pct = (len(matches) / total_targets) * 100
                    out.write(f"  - {primer_id} : {match_pct:.2f}% ({len(matches)}/{total_targets})\n")
                
                # Intersection brute (toutes les amorces présentes matchent)
                # Raw intersection (all present primers must match)
                if set_matches_list:
                    intersection_matches = set.intersection(*set_matches_list)
                else:
                    intersection_matches = set()
                    
                # Détermination des amorces essentielles présentes dans le set d'amorces
                # Determination of essential primers present in the primer set
                if args.pcr:
                    essential_in_set = [p for p in primers.keys() if p in ESSENTIAL_PCR]
                else:
                    essential_in_set = [p for p in primers.keys() if p in ESSENTIAL_LAMP]
                    
                # Fallback sur toutes les amorces s'il n'y a pas d'essentielles identifiées
                # Fallback to all primers if no essential primers are identified
                if not essential_in_set:
                    essential_in_set = list(primers.keys())
                    
                # Intersection essentielle (uniquement sur les amorces essentielles présentes)
                # Essential intersection (only on present essential primers)
                essential_matches_list = [primer_matches[set_id][p] for p in essential_in_set]
                if essential_matches_list:
                    essential_intersection_matches = set.intersection(*essential_matches_list)
                else:
                    essential_intersection_matches = set()
                    
                # Sélection de l'intersection pour la validation (brute ou essentielle)
                # Selection of the intersection for validation (raw or essential)
                if args.strict_intersection:
                    validation_matches = intersection_matches
                else:
                    validation_matches = essential_intersection_matches
                    
                valid_order_matches = []
                seq_details = []
                
                # Vérification de l'ordre et calcul de la taille de l'amplicon
                # Order verification and amplicon size calculation
                for seq_id in validation_matches:
                    positions = primer_positions[set_id][seq_id]
                    
                    # Tri des noms d'amorces selon la coordonnée de départ
                    sorted_primers = sorted(positions.keys(), key=lambda p: positions[p][0])
                    
                    expected_sense = [p for p in MASTER_ORDER if p in positions]
                    expected_anti = [p for p in MASTER_ORDER_RC if p in positions]
                    
                    is_correct_order = (sorted_primers == expected_sense) or (sorted_primers == expected_anti)
                    
                    if is_correct_order:
                        valid_order_matches.append(seq_id)
                        
                    starts = [pos[0] for pos in positions.values()]
                    ends = [pos[1] for pos in positions.values()]
                    amplicon_size = max(ends) - min(starts)
                    
                    status = txt['order_correct'] if is_correct_order else txt['order_incorrect']
                    observed_order = "-".join(sorted_primers)
                    seq_details.append((seq_id, amplicon_size, status, observed_order))
                    
                # Sauvegarde pour combine
                valid_sequences_per_set[set_id] = set(valid_order_matches)
                
                global_match_pct = (len(intersection_matches) / total_targets) * 100
                valid_match_pct = (len(valid_order_matches) / total_targets) * 100 if total_targets > 0 else 0
                
                out.write(f"\n{txt['global_raw']} : {global_match_pct:.2f}% ({len(intersection_matches)}/{total_targets})\n")
                out.write(f"{txt['global_valid']} : {valid_match_pct:.2f}% ({len(valid_order_matches)}/{total_targets})\n")
                
                # 1. Option : Ne pas afficher les séquences (si --summary-only)
                if not args.summary_only:
                    out.write(f"\n{txt['amplified_seqs'].format(set_id)}\n")
                    out.write(f"{txt['table_header']}\n")
                    
                    if seq_details:
                        seq_details.sort(key=lambda x: x[0])
                        for detail in seq_details:
                            out.write(f"{detail[0]}\t{detail[1]}\t{detail[2]}\t{detail[3]}\n")
                    else:
                        out.write(f"{txt['no_seq']}\n")
                
                # 2. Option : Exporter les séquences dans un fichier propre (--export-seqs)
                if args.export_seqs:
                    base_name, ext = os.path.splitext(args.output)
                    export_filename = f"{base_name}_Set{set_id}_seqs.txt"
                    try:
                        with open(export_filename, 'w') as f_export:
                            sorted_valid = sorted(valid_order_matches)
                            for s in sorted_valid:
                                f_export.write(f"{s}\n")
                        out.write("\n" + txt['export_success'].format(export_filename) + "\n")
                    except Exception as e:
                        out.write("\n" + txt['export_err'].format(e) + "\n")
                
                out.write("\n" + "="*40 + "\n\n")

            # 3. Option : Calcul Combinatoire (--combine)
            if args.combine and len(primer_sets) >= 2:
                out.write(f"{txt['combine_title']}\n")
                out.write("-" * len(txt['combine_title']) + "\n")
                out.write(f"{txt['combine_subtitle']}\n\n")
                
                set_keys = list(primer_sets.keys())
                for s1, s2 in itertools.combinations(set_keys, 2):
                    union_set = valid_sequences_per_set[s1] | valid_sequences_per_set[s2]
                    combine_pct = (len(union_set) / total_targets) * 100 if total_targets > 0 else 0
                    out.write(f"  - Set {s1} + Set {s2} : {combine_pct:.2f}% ({len(union_set)}/{total_targets})\n")
                out.write("\n" + "="*40 + "\n\n")

        print(txt['done'].format(args.output))
        
    except Exception as e:
        print(txt['write_err'].format(e))

if __name__ == "__main__":
    main()
