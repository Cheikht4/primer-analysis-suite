# Primer Aligner (Outil d'Alignement d'Amorces)

Ce script Python permet d'aligner automatiquement une liste d'amorces (primers) sur une séquence génomique cible ou sur un alignement multiple de séquences (MSA). Il génère en sortie un fichier FASTA où chaque amorce est parfaitement positionnée par rapport à la séquence cible, remplie de tirets (`-`) pour respecter l'alignement global.

Contrairement aux outils d'alignement multiple classiques (Clustal, MUSCLE) qui tentent d'aligner les amorces entre elles, cet outil traite chaque amorce individuellement par rapport à une séquence de référence.

## Fonctionnalités Principales

- **Support des génomes de toutes tailles** : Recherche efficace basée sur les expressions régulières.
- **Bases Dégénérées (IUPAC)** : Gère automatiquement les bases dégénérées (ex: `N`, `R`, `Y`, `W`, etc.) présentes dans vos amorces.
- **Tolérance aux erreurs** : Possibilité de tolérer un certain nombre de mésappariements (mismatches) ou d'indels entre l'amorce et la cible.
- **Support du brin antisens (Reverse Complement)** : L'outil cherche sur les deux brins. Si une amorce est trouvée sur le brin antisens, elle est automatiquement renommée avec le suffixe `+c` et sa séquence reverse-complémentaire est alignée.
- **Support des Alignements Multiples (MSA)** : Si vous fournissez un fichier contenant plusieurs séquences déjà alignées (avec des tirets `-`), l'outil utilisera la première séquence comme référence. Il injectera ensuite les mêmes trous (`-`) dans l'amorce pour ne pas casser les colonnes de votre alignement. La séquence d'origine de l'amorce n'est pas modifiée.

## Prérequis

- Python 3.6 ou supérieur.
- Bibliothèques Python : `biopython` et `regex`.

### Installation des dépendances

```bash
pip install biopython regex
```
*(ou `pip3 install biopython regex` selon votre système)*

## Utilisation

Le programme s'exécute en ligne de commande.

```bash
python3 align_primers.py -t <fichier_cible> -p <fichier_amorces> -o <fichier_sortie> [-e <erreurs_max>]
```

### Paramètres

- `-t` / `--target` **(Requis)** : Le fichier FASTA contenant votre génome cible OU votre alignement multiple (MSA). Dans le cas d'un MSA, la première séquence du fichier sera considérée comme la référence pour l'alignement des amorces.
- `-p` / `--primers` **(Requis)** : Le fichier contenant vos amorces. Il peut s'agir d'un fichier **FASTA** ou d'un simple fichier **Texte**.
  - Si vous utilisez un fichier texte, vous pouvez mettre une séquence par ligne, ou le format `Nom Sequence` (séparés par un espace ou une tabulation). S'il n'y a pas de nom, le script générera `Primer_1`, `Primer_2`, etc.
- `-o` / `--output` **(Requis)** : Le nom du fichier FASTA d'alignement qui sera généré.
- `-e` / `--errors` *(Optionnel)* : Le nombre maximal d'erreurs tolérées (mésappariements, insertions, délétions). La valeur par défaut est `2`. Mettez `0` pour exiger une correspondance exacte (aux bases dégénérées près).

## Exemples

### 1. Alignement sur une séquence simple
Vous avez un génome `genome.fasta` et des amorces `primers.fasta`.

```bash
python3 align_primers.py -t genome.fasta -p primers.fasta -o alignement_final.fasta
```
Le fichier `alignement_final.fasta` contiendra votre génome, suivi de vos amorces alignées en dessous avec des tirets.

### 2. Alignement sur un alignement multiple (MSA)
Vous avez un fichier `variants_alignes.fasta` contenant 10 séquences de virus alignées avec des gaps, et vous voulez voir où tombent vos amorces.

```bash
python3 align_primers.py -t variants_alignes.fasta -p primers.fasta -o variants_avec_amorces.fasta -e 1
```
Le fichier de sortie contiendra les 10 séquences d'origine parfaitement alignées, suivies des amorces. Si une amorce tombe sur une zone où l'alignement possède une insertion (gap dans la référence), l'outil insérera automatiquement un gap dans l'amorce au bon endroit.

## Remarques
- Les fichiers de sortie peuvent être ouverts dans n'importe quel visualiseur d'alignement standard (AliView, MEGA, Jalview, Geneious, etc.).
- Le programme ne modifie jamais les bases de vos amorces. S'il y a un mismatch entre l'amorce et le génome, la base de l'amorce sera affichée.
