# virtualPCR

**In silico PCR for simple and complex tasks**

[![Java](https://img.shields.io/badge/Java-25+-orange.svg)](https://www.oracle.com/java/technologies/downloads/)
[![Platform](https://img.shields.io/badge/Platform-Windows%20%7C%20Linux%20%7C%20macOS-blue.svg)]()
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)
[![Online Tool](https://img.shields.io/badge/Try%20Online-virtualPCR-green.svg)](https://primerdigital.com/tools/epcr.html)
[![DOI](https://img.shields.io/badge/DOI-10.3389%2Ffbinf.2024.1464197-blue.svg)](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2024.1464197/full)

---

## Table of Contents

- [Overview](#overview)
- [Citation](#citation)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Configuration File](#configuration-file)
- [Option Reference](#option-reference)
- [Input Formats](#input-formats)
- [Output](#output)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Author & Contact](#author--contact)

---

## Overview

virtualPCR is a versatile tool for **in silico PCR analysis**, designed to ensure primer/probe specificity across a wide range of applications. It enables researchers to predict primer/probe binding sites, assess mismatch tolerance, evaluate DNA duplex stability, and perform genome-wide searches.

### Applications

- Gene discovery via homology analysis
- Molecular diagnostics and primer validation
- Genome profiling and repeat sequence identification
- CRISPR-Cas gRNA target evaluation
- miRNA and probe specificity testing

### Key Capabilities

- Genome-wide primer/probe specificity analysis
- Linear and circular DNA support (plasmids, mitochondria, plastids)
- Bisulfite-treated DNA simulation (C→T conversion for methylation studies)
- Linked search mode for complex primer arrangements
- Degenerate nucleotides (IUPAC), LNA, and inosine support
- Batch file processing and automation

---

## Citation

If you use virtualPCR in your work, please cite:

> Kalendar R, Shevtsov A, Otarbay Z, Ismailova A. (2024). *In silico PCR analysis: a comprehensive bioinformatics tool for enhancing nucleic acid amplification assays*. **Frontiers in Bioinformatics**, 4:1464197. [DOI:10.3389/fbinf.2024.1464197](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2024.1464197/full)

---

## Requirements

| Requirement | Details |
|-------------|---------|
| **Java** | Version 25 or higher |
| **OS** | Windows, Linux, or macOS |
| **RAM** | Default JVM heap is sufficient for most tasks; increase for whole-genome analyses (see [Large Genomes](#large-genomes)) |

- **Download Java:** <https://www.oracle.com/java/technologies/downloads/>
- **Set Java Path:** <https://www.java.com/en/download/help/path.html>

---

## Installation

### Option 1: Direct Download

1. Download `virtualPCR.jar` from the `dist` directory of this repository.
2. Place it in your preferred location.
3. Ensure Java 25+ is installed and available in your `PATH`.

Verify Java is available:

```bash
java -version
```

### Option 2: Install Java via Conda

```bash
# Add conda-forge channel and set priority
conda config --add channels conda-forge
conda config --set channel_priority strict

# Create environment with OpenJDK 26
conda create -n java26 openjdk=26

# Activate environment
conda activate java26

# Verify installation
java -version
```

### Genome Downloads

Download individual chromosome FASTA files from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/genome/). For example, the human genome (T2T-CHM13v2.0) is available at [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/).

---

## Quick Start

```bash
# Basic usage
java -jar virtualPCR.jar config.file
```

### Getting Help

Run with no arguments — or with a help flag — to print usage and the full list of configuration options:

```bash
# Show help (any of these are equivalent)
java -jar virtualPCR.jar -help
java -jar virtualPCR.jar /?
```

Recognized help flags (case-insensitive): `-help`, `--help`, `-h`, `/help`, `/h`, `/?`, `-?`, `?`. Running with no configuration file also prints this help.

### Large Genomes

For genome-wide analyses, allocate additional memory using JVM flags:

```bash
# Linux / macOS
java -Xms4g -Xmx16g -jar virtualPCR.jar config.file

# Windows
java -Xms4g -Xmx16g -jar C:\virtualPCR\dist\virtualPCR.jar C:\virtualPCR\test\config.file
```

- **`-Xms`** — initial heap size (memory allocated at startup)
- **`-Xmx`** — maximum heap size (upper memory limit)

Increase `-Xmx` if you see `OutOfMemoryError`. A good rule of thumb for vertebrate genomes is `-Xmx16g` or higher.

---

## Configuration File

All parameters are specified in a plain-text configuration file (e.g., `config.file`), one `key=value` per line. Boolean options accept `true` or `false`.

```ini
# --- Paths ---
targets_path=C:\virtualPCR\test\ch02.fasta
primers_path=C:\virtualPCR\test\rt.txt
output_path=

# --- Search mode ---
type=primer
molecular=linear
linkedsearch=false
FRpairs=false
CTconversion=false

# --- Filters ---
minlen=200
maxlen=500
number3errors=1

# --- Output control ---
primerstatistic=true
SequenceExtract=true
ShowPrimerAlignment=true
ShowPCRProducts=true
ShowOnlyAmplicons=false
ShowPrimerAlignmentPCRproduct=false
```

> **Note:** Option keys are case-insensitive (`linkedsearch` and `LinkedSearch` are equivalent). For clarity, this README uses the casing shown in the configuration example above.

### Paths

| Parameter | Description |
|-----------|-------------|
| `targets_path` | Path to the target sequence file (FASTA or plain text). |
| `primers_path` | Path to the primer/probe file (FASTA or tab/space-delimited). |
| `output_path` | Output directory. Leave empty to write results next to the input. |

---

## Option Reference

### Summary

| Option | Values | Default | Purpose |
|--------|--------|---------|---------|
| `type` | `primer`, `probe` | `primer` | Search mode |
| `molecular` | `linear`, `circle` | `linear` | Template topology |
| `linkedsearch` | `true`, `false` | `false` | Enable linked/associated search |
| `FRpairs` | `true`, `false` | `false` | Restrict analysis to defined F/R pairs |
| `CTconversion` | `true`, `false` | `false` | Simulate bisulfite conversion |
| `minlen` | integer (bp) | `30` | Minimum amplicon length (inclusive) |
| `maxlen` | integer (bp) | `3000` | Maximum amplicon length (inclusive) |
| `number3errors` | integer | `1` | Allowed mismatches near the 3′ end |
| `primerstatistic` | `true`, `false` | `true` | Print per-primer summary statistics |
| `SequenceExtract` | `true`, `false` | `false` | Extract amplicon sequences into output |
| `ShowPrimerAlignment` | `true`, `false` | `true` | Show all primer-to-target alignments |
| `ShowPCRProducts` | `true`, `false` | `true` | Report predicted PCR products |
| `ShowOnlyAmplicons` | `true`, `false` | `false` | Report amplicon lengths only (no alignments) |
| `ShowPrimerAlignmentPCRproduct` | `true`, `false` | `true` | Restrict alignments to those forming products |

### `type` — Search Mode

| Value | Description |
|-------|-------------|
| `primer` | Standard PCR primer search **(default)**. |
| `probe` | Search for binding sites of primers, probes, and short sequences — TaqMan, Molecular Beacons, miRNA, CRISPR-Cas gRNA, microarray oligos, etc. Recommended when primer binding sites are not found, or when complementarity is expected for only part of the sequence (e.g., Molecular Beacon stems). |

### `molecular` — Template Topology

| Value | Description |
|-------|-------------|
| `linear` | Linear DNA **(default)**. |
| `circle` | Circular DNA (plasmids, mitochondrial/plastid genomes). Primers may produce one or two amplicons spanning the origin. |

### `minlen` / `maxlen` — Amplicon Size Filter

Defines the expected PCR product size range (in bp). Amplicons outside this range are filtered out.

- **Default:** `minlen=30`, `maxlen=3000`. Both bounds are **inclusive**; the valid range is 20–50000 bp.
- **Example:** `minlen=200` and `maxlen=500` restricts output to 200–500 bp products.

### `number3errors` — 3′ Mismatch Tolerance

Maximum number of mismatches allowed within the 3′ region of each primer. Lower values give stricter (PCR-realistic) specificity; higher values broaden the search.

### `CTconversion` — Bisulfite Conversion Simulation

Simulates bisulfite conversion for methylation studies. Only cytosines **not** followed by guanine (non-CpG) are converted to thymine on both strands:

```
5'-aaCGaagtCCCCa-3'        5'-aaCGaagtTTTTa-3'
   |||||||||||||     →         ||||||:|::::|
3'-ttGCttCaggggt-5'        3'-ttGCttTaggggt-5'
```

### `FRpairs` — Restrict to Defined Primer Pairs

When `true`, only the explicitly defined forward/reverse primer pairs are analyzed (no cross-pairing). Each line in the primer file defines one pair:

```
1020	aggcctgtgatgctgatgat	cccaacaccaaagaggaaag
1021	gctctgacctattgtctgtctgtct	cagtctccacagcagcagag
1022	gtcctgctgacacacaccact	cgcaggagtagaagaaaga
```

A single forward primer can be paired with multiple reverse primers (and vice versa) by listing additional columns on the same line.

### `linkedsearch` — Linked/Associated Search

A programmable search mode that locates binding sites for linked sequences within a specified distance. Supports tasks ranging from conventional sequence matching to in silico PCR with approximate matching.

**Syntax:** the first sequence is followed by the expected distance `[min-max]` and the second sequence:

```
>RT+(QMDVK)_RT-(YVDDML)
CARATGGAYGTNAARAC[200-300]TAYGTNGAYGAYATG
```

> The header `RT+(QMDVK)_RT-(YVDDML)` references conserved reverse-transcriptase protein motifs (QMDVK and YVDDML) used to derive the degenerate nucleotide primers.

**Variants:**

```
# Fixed distance (no range)
CARATGGAYGTNAARAC[300]TAYGTNGAYGAYATG

# Force the second sequence to be reverse-complemented (@ prefix)
CARATGGAYGTNAARAC[300]@CATRTCRTCNACRTA
```

### Output-control Options

| Option | Effect when `true` |
|--------|--------------------|
| `primerstatistic` | Prints a per-primer table (ID, sequence, binding-site hit count) sorted by frequency, plus a consolidated cross-file table and a companion `*.stats.tsv` file — see [Output](#output). |
| `SequenceExtract` | Writes extracted amplicon sequences into the output file. |
| `ShowPrimerAlignment` | Displays all stable primer-to-target alignments, including binding sites that may not produce PCR products. Useful for examining site stability, orientation, and coordinates. |
| `ShowPCRProducts` | Outputs predicted PCR products. Set to `false` to report binding sites only, without amplicon prediction. |
| `ShowOnlyAmplicons` | Outputs only amplicon lengths without detailed alignment analysis. Recommended for genome-wide in silico PCR with highly abundant repeat-based markers (iPBS, IRAP, ISSR, RAPD). |
| `ShowPrimerAlignmentPCRproduct` | Restricts alignment output to only those primer binding sites that contribute to predicted PCR products. |

---

## Input Formats

### Target Sequences

Sequence files can be plain text, single-entry FASTA, or multi-entry FASTA. No specific file extension is required. Template length is not limited (other than by available memory).

**FASTA format** consists of a header line starting with `>` followed by one or more lines of sequence:

```fasta
>sequence_id Optional description
ATCGATCGATCGATCGATCG...
```

### Primer/Probe Files

Primers can be provided in **tab-delimited**, **space-delimited**, or **FASTA** format. In delimited formats the first column is the primer ID and all subsequent columns are sequences belonging to that ID (one or more).

**Two-column tab-delimited:**

```
ITS1	TCCGTAGGTGAACCTGCGG
ITS2	GCTGCGTTCTTCATCGATGC
ITS3	GCATCGATGAAGAACGCAGC
ITS4	TCCTCCGCTTATTGATATGC
```

**Multi-column tab-delimited** (e.g. defined F/R pairs for use with `FRpairs=true`):

```
1020	aggcctgtgatgctgatgat	cccaacaccaaagaggaaag
1021	gctctgacctattgtctgtctgtct	cagtctccacagcagcagag
1022	gtcctgctgacacacaccact	cgcaggagtagaagaaaga
```

**FASTA:**

```fasta
>ITS1
TCCGTAGGTGAACCTGCGG
>ITS2
GCTGCGTTCTTCATCGATGC
>ITS3
GCATCGATGAAGAACGCAGC
>ITS4
TCCTCCGCTTATTGATATGC
```

**Rules:**

- Primer names may contain any characters, including spaces. Names need not be unique.
- Sequences must not contain spaces or non-nucleotide characters.
- Sequences shorter than 12 nt are skipped.
- File format (FASTA, tab, or space) is auto-detected.

### Supported Nucleotide Codes

| Code | Bases |
|------|-------|
| A | Adenine |
| C | Cytosine |
| G | Guanine |
| T | Thymine |
| U | Uracil (read as T) |
| I | Inosine (read as G) |
| R | A / G |
| Y | C / T |
| M | A / C |
| K | G / T |
| S | G / C |
| W | A / T |
| B | C / G / T |
| D | A / G / T |
| H | A / C / T |
| V | A / G / C |
| N | A / G / C / T |

**LNA (Locked Nucleic Acid) codes:** `dA = E`, `dC = F`, `dG = J`, `dT = L`

---

## Output

Results are saved as **tab-delimited, UTF-8 plain-text files**, containing (depending on options enabled):

- Primer binding-site coordinates and orientations
- Primer-target alignment details and mismatch analysis
- Predicted PCR product sizes
- Extracted amplicon sequences
- Per-primer summary statistics (when `primerstatistic=true`)

When a directory of target files is processed, all per-file reports are written to a single combined output file, followed by the global statistics described below.

### Per-primer statistics

With `primerstatistic=true`, each report ends with a per-primer table sorted by hit count (`PrimerID`, `Sequence`, `Hits`), followed by a **consolidated global table** that sums each primer's occurrence across every analyzed sequence and file:

| Column | Meaning |
|--------|---------|
| `Hits` | Total binding sites found |
| `Sequences` | Number of sequences in which the primer occurs |
| `Files` | Number of files in which the primer occurs |

When a directory of targets is processed, this table aggregates the whole batch. It is written both at the end of the main report and to a separate companion **`<report>.stats.tsv`** file (tab-separated, UTF-8) for direct import into spreadsheets or R.

---

## Examples

### 1. Standard primer screen against a chromosome

```ini
targets_path=test/ch02.fasta
primers_path=test/its.txt
type=primer
molecular=linear
minlen=100
maxlen=2000
number3errors=1
ShowPrimerAlignment=true
ShowPCRProducts=true
```

### 2. Genome-wide repeat-marker scan (iPBS / IRAP / ISSR)

```ini
targets_path=genome/chm13v2.fasta
primers_path=primers/ipbs.txt
type=primer
maxlen=5000
ShowOnlyAmplicons=true
```

Run with extra memory:

```bash
java -Xms4g -Xmx32g -jar virtualPCR.jar config.file
```

### 3. Probe / gRNA specificity check

```ini
targets_path=genome/hg38.fasta
primers_path=probes/grna.fasta
type=probe
ShowPCRProducts=false
ShowPrimerAlignment=true
```

### 4. Bisulfite-converted template

```ini
targets_path=test/bs_region.fasta
primers_path=test/bs_primers.txt
type=primer
CTconversion=true
minlen=80
maxlen=400
```

---

## Troubleshooting

| Symptom | Likely cause / fix |
|---------|-------------------|
| `UnsupportedClassVersionError` | Java version is older than 25. Install Java 25+ and verify with `java -version`. |
| `OutOfMemoryError` | Increase JVM heap with `-Xmx`, e.g. `-Xmx16g`. |
| No products reported | Try `type=probe`, increase `number3errors`, or widen `minlen`/`maxlen`. |
| Too many products on a genome | Tighten `number3errors`, narrow the amplicon size window, or switch to `ShowOnlyAmplicons=true` for compact output. |
| Primer file appears empty | Ensure sequences are ≥ 12 nt; check that ID and sequence are separated by a tab or space. |
| Not sure which options exist | Run `java -jar virtualPCR.jar -help` (or with no arguments) to print usage and all configuration options. |

---

## License

Released under the **GNU General Public License v3.0**. See [LICENSE.txt](LICENSE.txt) for the full text.

---

## Author & Contact

**Ruslan Kalendar**
📧 <ruslan.kalendar@helsinki.fi>

**Online version:** <https://primerdigital.com/tools/epcr.html>
