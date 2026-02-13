# virtualPCR

**In silico PCR for simple and complex tasks**

[![Java](https://img.shields.io/badge/Java-24+-orange.svg)](https://www.oracle.com/java/technologies/downloads/)
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
- [Author & Contact](#author--contact)

---

## Overview

virtualPCR is a versatile tool for conducting **in silico PCR analysis**, designed to ensure primer/probe specificity across a wide range of applications. It enables researchers to predict primer/probe binding sites, assess mismatch tolerance, evaluate DNA duplex stability, and perform genome-wide searches.

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

Kalendar R, Shevtsov A, Otarbay Z, Ismailova A. (2024). *In silico PCR analysis: a comprehensive bioinformatics tool for enhancing nucleic acid amplification assays*. **Frontiers in Bioinformatics**, 4:1464197.
[DOI:10.3389/fbinf.2024.1464197](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2024.1464197/full)

---

## Requirements

| Requirement | Details |
|-------------|---------|
| **Java** | Version 24 or higher |
| **OS** | Windows, Linux, or macOS |
| **RAM** | Default is sufficient; increase for large genomes (see [Memory](#large-genomes)) |

**Download Java:** https://www.oracle.com/java/technologies/downloads/
**Set Java Path:** https://www.java.com/en/download/help/path.html

---

## Installation

### Option 1: Direct Download

1. Download `virtualPCR.jar` from the `dist` directory
2. Place it in your preferred location
3. Ensure Java 24+ is installed and available in your PATH

### Option 2: Install Java via Conda

```bash
# Add conda-forge channel and set priority
conda config --add channels conda-forge
conda config --set channel_priority strict

# Create environment with OpenJDK 25
conda create -n java25 openjdk=25

# Activate environment
conda activate java25

# Verify installation
java -version
```

### Genome Downloads

Download individual chromosome FASTA files from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/genome/). For example, the human genome (T2T-CHM13v2.0) is available at [NCBI FTP](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/).

---

## Quick Start

```bash
# Basic usage
java -jar virtualPCR.jar config.file

# With increased memory for large genomes
java -Xms4g -Xmx16g -jar virtualPCR.jar config.file
```

### Large Genomes

For genome-wide analyses, allocate additional memory using JVM flags:

```bash
java -Xms4g -Xmx16g -jar C:\virtualPCR\dist\virtualPCR.jar C:\virtualPCR\test\config.file
```

- **`-Xms`** — initial heap size (memory allocated at startup)
- **`-Xmx`** — maximum heap size (upper memory limit)

---

## Configuration File

All parameters are specified in a plain text file (e.g., `config.file`):

```ini
targets_path=C:\virtualPCR\test\ch02.fasta
output_path=
primers_path=C:\virtualPCR\test\rt.txt
type=primer/probe
linkedsearch=false/true
molecular=linear/circle
primerstatistic=true
number3errors=1
minlen=200
maxlen=500
FRpairs=false/true
CTconversion=false/true
SequenceExtract=true
ShowPrimerAlignment=true
ShowOnlyAmplicons=true/false
ShowPCRProducts=true/false
ShowPrimerAlignmentPCRproduct=true/false
```

### Paths

| Parameter | Description |
|-----------|-------------|
| `targets_path` | Path to the target sequence file (FASTA or plain text) |
| `primers_path` | Path to the primer/probe file (FASTA or tab-delimited) |
| `output_path` | Output directory (leave empty for same directory as input) |

---

## Option Reference

### `type` — Search Mode

| Value | Description |
|-------|-------------|
| `primer` | Standard PCR primer search **(default)** |
| `probe` | Search for binding sites of primers, probes, and short sequences — TaqMan, Molecular Beacons, miRNA, CRISPR-Cas gRNA, microarray oligos, etc. Recommended when primer binding sites are not found, or when complementarity is expected for only part of the sequence (e.g., Molecular Beacon stems) |

### `molecular` — Template Topology

| Value | Description |
|-------|-------------|
| `linear` | Linear DNA **(default)** |
| `circle` | Circular DNA (plasmids, mitochondrial/plastid genomes). Primers may produce one or two amplicons spanning the origin |

### `minlen` / `maxlen` — Amplicon Size Filter

Defines the expected PCR product size range (in bp). Amplicons outside this range are filtered out.

- **Default:** 5000 bp maximum
- **Example:** `minlen=200` and `maxlen=500` restricts output to 200–500 bp products

### `CTconversion` — Bisulfite Conversion Simulation

Simulates bisulfite conversion for methylation studies. Only cytosines **not** followed by guanine (non-CpG) are converted to thymine on both strands:

```
5'-aaCGaagtCCCCa-3'        5'-aaCGaagtTTTTa-3'
   |||||||||||||     →         ||||||:|::::|
3'-ttGCttCaggggt-5'        3'-ttGCttTaggggt-5'
```

### `FRpairs` — Restrict to Defined Primer Pairs

When enabled, only defined forward/reverse primer pairs are analyzed. Each line in the primer file defines one pair (tab-separated):

```
1020	aggcctgtgatgctgatgat	cccaacaccaaagaggaaag
1021	gctctgacctattgtctgtctgtct	cagtctccacagcagcagag
1022	gtcctgctgacacacaccact	cgcaggagtagaagaaaga
```

A single forward primer can be paired with multiple reverse primers (and vice versa).

### `ShowPrimerAlignment` — Display All Binding Sites

When `true` (default), displays all stable primer-to-target alignments, including sites that may not produce PCR products under current conditions. Useful for examining binding site stability, orientation, and coordinates.

### `ShowPrimerAlignmentPCRproduct` — Alignments for PCR Products Only

When `true`, restricts alignment output to only those primer binding sites that contribute to predicted PCR products.

### `ShowPCRProducts` — PCR Product Prediction

When `true` (default), outputs predicted PCR products. Set to `false` to search for primer binding sites only, without amplicon prediction.

### `ShowOnlyAmplicons` — Amplicon Lengths Only

When `true`, outputs only amplicon lengths without detailed alignment analysis. Recommended for genome-wide in silico PCR with highly abundant repeat-based markers (iPBS, IRAP, ISSR, RAPD).

### `LinkedSearch` — Linked/Associated Search

A programmable search mode where binding sites are found for linked sequences within a specified distance. Linked searching supports tasks ranging from conventional sequence matching to in silico PCR with approximate matching.

**Syntax:** The forward primer sequence is followed by the expected distance `[min-max]` and the second sequence:

```
>RT+(QMDVK)_RT-(YVDDML)
CARATGGAYGTNAARAC[200-300]TAYGTNGAYGAYATG
```

**Variants:**

```
# Fixed distance (no range)
CARATGGAYGTNAARAC[300]TAYGTNGAYGAYATG

# Force second sequence to be reverse-complemented (@ prefix)
CARATGGAYGTNAARAC[300]@CATRTCRTCNACRTA
```

### `number3errors` — 3' Mismatch Tolerance

Number of allowed mismatches at the 3' end of the primer.

---

## Input Formats

### Target Sequences

Sequence files can be plain text, FASTA, or multi-entry FASTA. No file extension is required. Template length is not limited.

**FASTA format** consists of a header line starting with `>` followed by one or more lines of sequence:

```fasta
>sequence_id Optional description
ATCGATCGATCGATCGATCG...
```

### Primer/Probe Files

Primers can be provided in **tab-delimited** or **FASTA** format. The software reads the first two columns (name and sequence) or three/four columns for primer-probe sets.

**Tab-delimited format:**

```
ITS1	TCCGTAGGTGAACCTGCGG
ITS2	GCTGCGTTCTTCATCGATGC
ITS3	GCATCGATGAAGAACGCAGC
ITS4	TCCTCCGCTTATTGATATGC
```

**FASTA format:**

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

Primer names can contain any characters (including spaces). Names may be identical across entries. Sequences must not contain spaces or non-nucleotide characters.

### Supported Nucleotide Codes

| Code | Bases | Code | Bases |
|------|-------|------|-------|
| A | Adenine | M | A/C |
| T | Thymine | R | A/G |
| G | Guanine | W | A/T |
| C | Cytosine | S | G/C |
| U | Uracil | Y | C/T |
| I | Inosine | K | G/T |
| N | A/G/C/T | V | A/G/C |
| | | H | A/C/T |
| | | D | A/G/T |
| | | B | C/G/T |

**LNA (Locked Nucleic Acid) codes:** dA=E, dC=F, dG=J, dT=L

---

## Output

Results are saved as **tab-delimited plain text files**, containing (depending on options enabled):

- Primer binding site coordinates and orientations
- Primer-target alignment details and mismatch analysis
- Predicted PCR product sizes
- Extracted amplicon sequences

---

## Author & Contact

**Ruslan Kalendar**
📧 ruslan.kalendar@helsinki.fi

**Online Version:** https://primerdigital.com/tools/epcr.html
