#!/usr/bin/env python3
# =============================================================================
# primer_analyse.py - Pipeline unifié d'analyse d'amorces
# =============================================================================
# Auteur : Antigravity
# Description (FR) : Exécute align_primers.py, analyze_dimers_v2.py et lamp_coverage.py
# Description (EN) : Runs align_primers.py, analyze_dimers_v2.py and lamp_coverage.py
# =============================================================================

import argparse
import subprocess
import os
import sys
from Bio import SeqIO

# S'assurer que les imports locaux fonctionnent
# Ensure local imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def main():
    parser = argparse.ArgumentParser(description="Pipeline d'analyse d'amorces (Alignement, Dimères, Couverture) / Primer analysis pipeline")
    
    # Options générales
    parser.add_argument("-t", "--target", required=False, help="Fichier FASTA cible / Target FASTA file (Requis pour align et lamp)")
    parser.add_argument("-p", "--primers", required=True, help="Fichier FASTA ou TXT des amorces / Primers FASTA or TXT file")
    parser.add_argument("-o", "--output_prefix", default="results", help="Préfixe pour les fichiers de sortie / Prefix for output files")
    parser.add_argument("-e", "--errors", default="2", help="Nombre max d'erreurs tolérées / Max errors tolerated (défaut: 2)")
    parser.add_argument("--lng", choices=["en", "fr"], default="fr", help="Langue des rapports / Report language (défaut: fr)")
    parser.add_argument("--steps", nargs="+", choices=["align", "dimers", "lamp"], default=["align", "dimers", "lamp"], help="Analyses à exécuter / Analyses to run (défaut: toutes / all)")
    
    # Options spécifiques à lamp_coverage.py
    lamp_group = parser.add_argument_group("Options spécifiques à l'analyse de couverture (lamp_coverage.py)")
    lamp_group.add_argument("-s", "--strict-3prime", default="3", help="Taille zone 3' stricte / Strict 3' region size. Def: 3")
    lamp_group.add_argument("--summary-only", action="store_true", help="N'affiche que les statistiques (lamp) / Output only summary statistics")
    lamp_group.add_argument("--combine", action="store_true", help="Couverture combinatoire 2 à 2 (lamp) / Calculate 2-by-2 multiplexing coverage")
    lamp_group.add_argument("--strict-3prime-tolerate", type=int, choices=[0, 1, 2], default=0, help="Niveau de tolérance en zone 3' (0: tout strict, 1: pos 2 tolérée, 2: pos 1 et 2 tolérées) / Tolerance level in the 3' region (0, 1, or 2)")
    lamp_group.add_argument("--export-seqs", action="store_true", help="Exporte les séquences validées (lamp) / Export validated sequences per set")
    lamp_group.add_argument("--pcr", action="store_true", help="Mode PCR (lamp) / PCR mode")
    lamp_group.add_argument("--strict-intersection", action="store_true", help="Exige que toutes les amorces du fichier matchent la cible (lamp) / Requires all primers in the file to match the target")
    lamp_group.add_argument("--max-n-run", type=int, default=10, dest="max_n_run",
        help="Exclure les séquences avec un run de N consécutifs ≥ cette valeur / Exclude sequences with N-run >= this value. 0=désactivé/disabled. Def: 10")

    args = parser.parse_args()
    
    if ("align" in args.steps or "lamp" in args.steps) and not args.target:
        print("❌ Erreur : L'argument -t/--target est requis pour les analyses 'align' et 'lamp'.")
        sys.exit(1)
        
    if args.target and not os.path.exists(args.target):
        print(f"❌ Erreur : Le fichier cible '{args.target}' est introuvable.")
        sys.exit(1)
        
    if not os.path.exists(args.primers):
        print(f"❌ Erreur : Le fichier d'amorces '{args.primers}' est introuvable.")
        sys.exit(1)
        
    print(f"🚀 Début du pipeline 'primer_analyse'...")
    if args.target:
        print(f"Cible : {args.target}")
    print(f"Amorces : {args.primers}")
    print(f"Étapes actives : {', '.join(args.steps)}\n")
    
    aligned_fasta = None
    primers_fasta = None
    report_file = None
    
    # 1. Alignement des amorces / Align primers
    if "align" in args.steps:
        aligned_fasta = f"{args.output_prefix}_aligned.fasta"
        print("==================================================")
        print(f"🛠️  Étape 1 : Alignement des amorces (align_primers.py)")
        print("==================================================")
        cmd_align = [
            sys.executable, "align_primers.py",
            "-t", args.target,
            "-p", args.primers,
            "-o", aligned_fasta,
            "-e", str(args.errors)
        ]
        subprocess.run(cmd_align)
        print()
    
    # 2. Analyse des dimères / Dimer analysis
    if "dimers" in args.steps:
        print("==================================================")
        print(f"🧬 Étape 2 : Analyse des dimères (analyze_dimers_v2.py)")
        print("==================================================")
        
        primers_fasta = f"{args.output_prefix}_primers.fasta"
        try:
            import align_primers
            primers_records = align_primers.load_primers(args.primers)
            SeqIO.write(primers_records, primers_fasta, "fasta")
            
            input_data = f"{args.lng}\n{primers_fasta}\ny\n"
            cmd_dimers = [sys.executable, "analyze_dimers_v2.py"]
            subprocess.run(cmd_dimers, input=input_data, text=True)
            
        except Exception as e:
            print(f"⚠️ Erreur lors de l'exécution de l'analyse des dimères / Error during dimer analysis: {e}")
        print()
        
    # 3. Couverture LAMP / LAMP Coverage
    if "lamp" in args.steps:
        print("==================================================")
        print(f"🔍 Étape 3 : Couverture LAMP (lamp_coverage.py)")
        print("==================================================")
        report_file = f"{args.output_prefix}_coverage_report.txt"
        
        cmd_lamp = [
            sys.executable, "lamp_coverage.py",
            "-t", args.target,
            "-p", args.primers,
            "-o", report_file,
            "-e", str(args.errors),
            "-s", str(args.strict_3prime),
            "--lng", args.lng
        ]
        
        if args.summary_only:
            cmd_lamp.append("--summary-only")
        if args.combine:
            cmd_lamp.append("--combine")
        if args.export_seqs:
            cmd_lamp.append("--export-seqs")
        if args.pcr:
            cmd_lamp.append("--pcr")
        if args.strict_3prime_tolerate > 0:
            cmd_lamp.extend(["--strict-3prime-tolerate", str(args.strict_3prime_tolerate)])
        if args.strict_intersection:
            cmd_lamp.append("--strict-intersection")
        if args.max_n_run != 10:  # Ne transmet que si différent de la valeur par défaut / Only pass if different from default
            cmd_lamp.extend(["--max-n-run", str(args.max_n_run)])
            
        subprocess.run(cmd_lamp)
        print()
    
    print(f"✅ Pipeline terminé avec succès ! / Pipeline completed successfully!")
    print(f"Fichiers générés / Generated files:")
    if primers_fasta:
        print(f"  - FASTA des amorces converties / Converted primers FASTA : {primers_fasta}")
    if "align" in args.steps and aligned_fasta:
        print(f"  - Alignement MSA / MSA Alignment : {aligned_fasta}")
    if "dimers" in args.steps:
        print(f"  - Rapports dimères / Dimer reports : dimer_analysis_{args.output_prefix}_primers (.txt, .csv, .html)")
    if "lamp" in args.steps and report_file:
        print(f"  - Couverture LAMP / LAMP Coverage : {report_file}")

if __name__ == "__main__":
    main()
