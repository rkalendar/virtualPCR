# virtualPCR

## In silico PCR for simple and complex tasks
**Author:** Ruslan Kalendar  
**Email:** ruslan.kalendar@helsinki.fi

virtualPCR is a versatile software tool for conducting **in silico PCR analysis**, designed to ensure primer/probe specificity across a wide range of applications. It enables researchers to predict primer/probe binding sites, assess mismatch tolerance, evaluate DNA duplex stability, and perform genome-wide searches. This makes it valuable for:

- Gene discovery via homology analysis
- Molecular diagnostics
- Genome profiling
- Repeat sequence identification
- CRISPR-Cas gRNA target evaluation
- miRNA and probe specificity testing

The tool supports **linear and circular DNA templates**, including **bisulfite-treated DNA**, and allows batch file processing for large datasets.

---

## Citation
Kalendar R, Shevtsov A, Otarbay Z, Ismailova A. (2024). *In silico PCR analysis: a comprehensive bioinformatics tool for enhancing nucleic acid amplification assays*. **Frontiers in Bioinformatics**, 4:1464197. [DOI:10.3389/fbinf.2024.1464197](https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2024.1464197/full)

---

## Online Access
Use virtualPCR online at: [https://primerdigital.com/tools/epcr.html](https://primerdigital.com/tools/epcr.html)

---

## Installation & Requirements
- **Programming Language:** Java 24+
- **Java Downloads:** [Oracle Java](https://www.oracle.com/java/technologies/downloads/)
- **Set Java Path:** [Instructions](https://www.java.com/en/download/help/path.html)

### Genome Downloads
If analyzing a specific genome, download individual chromosome FASTA files from [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/genome/). For example, the **human genome (T2T-CHM13v2.0)** is available at: [NCBI FTP](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/).

---

## Running virtualPCR
The main executable is `virtualPCR.jar`, located in the `dist` directory. Copy it to any location for use.

### Example Command
```bash
java -jar C:\virtualPCR\dist\virtualPCR.jar C:\virtualPCR\test\config.file
```

### Large Genomes (Increase Memory)
```bash
java -Xms4g -Xmx16g -jar C:\virtualPCR\dist\virtualPCR.jar C:\virtualPCR\test\config.file
```

---

## Configuration File
Parameters are specified in a plain text file (e.g., `config.file`):
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

### Key Options
- **type=probe** → search for primers and probes (TaqMan, Molecular Beacon, miRNA, CRISPR, microarrays)
- **minlen / maxlen** → define expected PCR product size (default: 5000 bp)
- **ShowPCRProducts** → enable/disable predicted PCR product output
- **molecular** → analyze linear or circular DNA (plasmids, mitochondria, plastids)
- **CTconversion** → simulate bisulfite conversion for methylation studies
- **FRpairs** → restrict search to defined forward/reverse primer pairs
- **ShowPrimerAlignment** → show all stable primer-target alignments
- **ShowOnlyAmplicons** → output only amplicon sizes (recommended for genome-wide repeat-based PCR)
- **LinkedSearch** → search linked primer sites with defined distance ranges

---

## Input & Output
- **Input:** Sequences in plain text, FASTA, or tab-delimited format. Supports degenerate nucleotides (IUPAC codes) and special cases (LNA, inosine).
- **Output:** Results saved as tab-delimited plain text files.

### Example Primer File (Tabular)
```tsv
ITS1    TCCGTAGGTGAACCTGCGG
ITS2    GCTGCGTTCTTCATCGATGC
ITS3    GCATCGATGAAGAACGCAGC
ITS4    TCCTCCGCTTATTGATATGC
```

### Example Primer File (FASTA)
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

---

## Features at a Glance
- 🔹 Genome-wide primer/probe specificity analysis
- 🔹 Linear & circular DNA support
- 🔹 Bisulfite DNA simulation (C>T conversion)
- 🔹 Linked search mode for complex primer arrangements
- 🔹 Batch file processing & automation
- 🔹 User-friendly configuration via text files

---

## License & Contact
For academic and research use.  
Contact: **ruslan.kalendar@helsinki.fi**

