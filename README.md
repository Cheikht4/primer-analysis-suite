# 🧬 Primer Analysis Suite

A comprehensive bioinformatics pipeline for **LAMP and PCR primer analysis**, including multiple sequence alignment (MSA), primer-dimer detection, and genome coverage evaluation.

> Developed for Dengue virus primer validation, but applicable to any viral or bacterial target genome.

---

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Project Structure](#project-structure)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
  - [Unified Pipeline (`primer_analyse.py`)](#unified-pipeline)
  - [Alignment only (`align_primers.py`)](#alignment-only)
  - [Dimer Analysis (`analyze_dimers_v2.py`)](#dimer-analysis)
  - [LAMP Coverage (`lamp_coverage.py`)](#lamp-coverage)
- [Examples](#examples)
- [Output Files](#output-files)
- [Notes](#notes)

---

## 🔍 Overview

This suite of Python scripts provides a modular, command-line pipeline to:

1. **Align** primer sequences onto a reference genome or a Multiple Sequence Alignment (MSA)
2. **Detect** primer-dimer interactions (hairpin, homodimer, heterodimer) using `ntthal` (Primer3)
3. **Evaluate** the coverage of your primer sets across a population of viral genome sequences

---

## ✨ Features

- ✅ **Supports all genome sizes** — efficient regex-based search
- ✅ **Degenerate/IUPAC bases** — handles `N`, `R`, `Y`, `W`, `S`, `M`, `K`, `B`, `D`, `H`, `V` automatically
- ✅ **Error tolerance** — configurable number of mismatches/indels
- ✅ **Reverse complement search** — searches both strands; antisense primers are auto-detected and renamed with `+c`
- ✅ **MSA-compatible** — injects gaps into primers to preserve alignment columns
- ✅ **Dimer detection** — hairpin, homodimer, and heterodimer analysis via `ntthal`
- ✅ **LAMP coverage** — evaluates per-sequence coverage with strict 3' zone validation
- ✅ **Combinatorial mode** — tests all pairwise combinations of primer sets for coverage
- ✅ **Multilingual** — reports available in French (`fr`) and English (`en`)
- ✅ **Modular pipeline** — run all steps or only the ones you need with `--steps`

---

## 📁 Project Structure

```
.
├── primer_analyse.py       # 🚀 Main unified pipeline (entry point)
├── align_primers.py        # Step 1: MSA primer alignment
├── analyze_dimers_v2.py    # Step 2: Primer-dimer analysis
├── lamp_coverage.py        # Step 3: LAMP/PCR genome coverage
├── requirements.txt        # Python dependencies
└── README.md               # This file
```

---

## 🛠️ Prerequisites

- **Python 3.6+**
- **Primer3** (for dimer analysis — provides `ntthal`)
  - macOS: `brew install primer3`
  - Linux (Ubuntu/Debian): `sudo apt-get install primer3`
- **Python libraries**: `biopython`, `regex`, `tqdm`

---

## 💻 Installation

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/primer-analysis-suite.git
cd primer-analysis-suite
```

### 2. Install Python dependencies

```bash
pip install -r requirements.txt
```

### 3. Install Primer3 (required for dimer analysis)

```bash
# macOS
brew install primer3

# Linux (Ubuntu/Debian)
sudo apt-get install primer3
```

Verify installation:
```bash
ntthal -h
```

---

## 🚀 Usage

### Unified Pipeline

`primer_analyse.py` is the **main entry point**. It orchestrates all three analysis steps.

```bash
python3 primer_analyse.py -t <target.fasta> -p <primers.txt> [options]
```

#### General Arguments

| Argument | Short | Description | Default |
|---|---|---|---|
| `--target` | `-t` | Target FASTA file (genome or MSA) | *Required* |
| `--primers` | `-p` | Primers file (FASTA or TXT format) | *Required* |
| `--output_prefix` | `-o` | Prefix for all output files | `results` |
| `--errors` | `-e` | Max allowed mismatches/indels | `2` |
| `--lng` | | Report language (`fr` or `en`) | `fr` |
| `--steps` | | Steps to run: `align`, `dimers`, `lamp` | all three |

#### LAMP Coverage Options (`--steps lamp`)

| Argument | Short | Description | Default |
|---|---|---|---|
| `--strict-3prime` | `-s` | Size of the strict 3' region (no errors allowed) | `3` |
| `--summary-only` | | Output only summary statistics | off |
| `--combine` | | Test all pairwise primer set combinations | off |
| `--export-seqs` | | Export validated sequences per primer set | off |
| `--pcr` | | Switch to PCR mode (instead of LAMP) | off |

---

### Alignment Only

Run primer alignment independently:

```bash
python3 align_primers.py -t <target.fasta> -p <primers.txt> -o <output.fasta> [-e <max_errors>]
```

**Arguments:**
- `-t` / `--target` : FASTA file of the target genome or MSA *(required)*
- `-p` / `--primers` : Primers file (FASTA or plain text, one per line) *(required)*
- `-o` / `--output` : Output FASTA alignment file *(required)*
- `-e` / `--errors` : Max errors tolerated (default: `2`)

**Primers file format (TXT):**
```
# One sequence per line, or "Name Sequence" format
F3  ATGCATGCATGC
B3  CGTACGTACGTA
FIP ATGCATGCATGCATGCATGC
BIP CGTACGTACGTACGTACGTA
```

---

### Dimer Analysis

Interactive dimer analysis (can also be invoked via the unified pipeline):

```bash
python3 analyze_dimers_v2.py
```

The script will prompt for:
1. Language (`fr` or `en`)
2. Path to your primers FASTA file
3. Whether to include all-vs-all heterodimer analysis

**Detects:**
- 🔁 **Hairpin** structures (self-folding)
- 🔗 **Homodimers** (primer with itself)
- ⚡ **Heterodimers** (primer pairs)

Results are saved as `.txt`, `.csv`, and `.html` reports.

---

### LAMP Coverage

Evaluate how well your primer set covers a population of genome sequences:

```bash
python3 lamp_coverage.py -t <sequences.fasta> -p <primers.txt> -o <report.txt> [options]
```

**Arguments:**
- `-t` / `--target` : FASTA file with all sequences to test
- `-p` / `--primers` : Primers file
- `-o` / `--output` : Output report file
- `-e` / `--errors` : Max mismatches outside strict 3' zone (default: `2`)
- `-s` / `--strict-3prime` : Size of the strict 3' region (no errors, default: `3`)
- `--combine` : Test pairwise combinations of primer sets
- `--summary-only` : Only output summary stats (no per-sequence detail)
- `--export-seqs` : Export sequences passing each primer set to separate files
- `--pcr` : Use PCR mode
- `--lng` : Language of report (`fr` or `en`)

---

## 📊 Examples

### Run the full pipeline

```bash
python3 primer_analyse.py \
    -t sequences_DENGUE_1.fasta \
    -p primer_dengue_1.txt \
    -o dengue1_results \
    -e 2 \
    --lng fr
```

### Run only LAMP coverage with combinatorial mode

```bash
python3 primer_analyse.py \
    -t sequences_DENGUE_1.fasta \
    -p primer_dengue_1.txt \
    -o dengue1_results \
    --steps lamp \
    --combine \
    --summary-only \
    -s 5
```

### Run alignment and dimer analysis only (no coverage)

```bash
python3 primer_analyse.py \
    -t sequences_DENGUE_1.fasta \
    -p primer_dengue_1.txt \
    -o dengue1_results \
    --steps align dimers
```

### Strict exact alignment (no errors)

```bash
python3 align_primers.py \
    -t genome.fasta \
    -p primers.fasta \
    -o strict_alignment.fasta \
    -e 0
```

---

## 📂 Output Files

| File | Generated by | Description |
|---|---|---|
| `<prefix>_aligned.fasta` | `align` step | MSA FASTA with primers inserted |
| `<prefix>_primers.fasta` | `dimers` step | Primers converted to FASTA |
| `dimer_analysis_<prefix>_primers.txt` | `dimers` step | Full dimer analysis report |
| `dimer_analysis_<prefix>_primers.csv` | `dimers` step | Machine-readable dimer data |
| `dimer_analysis_<prefix>_primers.html` | `dimers` step | Interactive HTML dimer report |
| `<prefix>_coverage_report.txt` | `lamp` step | Per-sequence coverage report |

---

## 📝 Notes

- Output alignment files can be opened in any standard alignment viewer: **AliView**, **MEGA**, **Jalview**, **Geneious**, etc.
- The scripts **never modify your primer bases**. If there is a mismatch between a primer and the genome, the primer's original base is shown.
- Antisense primers found on the reverse strand are automatically renamed with the `+c` suffix and their reverse-complement sequence is aligned.
- For large FASTA files (> 1000 sequences), use `--summary-only` to avoid very large report files.

---

## 📄 License

This project is released under the MIT License.

---

## 👤 Author

**Cheikh Talibou Ya** — Primer analysis pipeline for viral genome surveillance.
