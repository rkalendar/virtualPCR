## virtualPCR 
## "In silico PCR for simple and complex tasks"
by Ruslan Kalendar 
Email: ruslan.kalendar@helsinki.fi

In silico PCR analysis is a valuable and productive adjunctive approach for ensuring primer or probe specificity across a broad spectrum of PCR applications encompassing gene discovery through homology analysis, molecular diagnostics, DNA profiling, and repeat sequence identification. The prediction of primer, probe and microRNA (miRNA) sensitivity and specificity necessitates thorough database searches, accounting for an optimal balance of mismatch tolerance, sequence similarity, and thermal stability. This software facilitates in silico PCR analyses of both linear and circular DNA templates, including bisulfited treatment DNA, enabling multiple primer or probe searches within databases of varying scales alongside advanced search functionalities. This tool is suitable for processing batch files and is essential for automation when working with large amounts of data.

## Citation reference
Kalendar R, Shevtsov A, Otarbay Z, Ismailova A 2024. In silico PCR analysis: a comprehensive bioinformatics tool for enhancing nucleic acid amplification assays. Frontiers in Bioinformatics, 4: 1464197. DOI:10.3389/fbinf.2024.1464197
https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2024.1464197/full

## Online virtualPCR: 
https://primerdigital.com/tools/epcr.html

## Availability and requirements:
Programming language: Java 23 or higher
Java Downloads: 
https://www.oracle.com/java/technologies/downloads/

How do I set or change the Java path system variable: 
https://www.java.com/en/download/help/path.html

To run the project from the command line. Command-line options, separated by spaces. 
The executive file ```virtualPCR.jar``` is in the ```dist``` directory, which can be copied to any location. 
The program needs to specify a file containing all the necessary parameters and paths to the files to be analysed. Here for an example the text file ```config.file``` is given. It is necessary to specify the path to this file:

```java -jar virtualPCR.jar <config.file_path/config.file>```

```java -jar C:\virtualPCR\dist\virtualPCR.jar C:\virtualPCR\test\config.file```


### Basic usage: 
To enter parameters and specify the location of the target files and primer's file, you must specify this via a file on the command line. An example of such a file here ```config.file``` (file name or extension does not matter):

> **config.file**
```
targets_path=C:\virtualPCR\test\1.txt
primers_path=C:\virtualPCR\task\primers.txt
type=primer/probe
linkedsearch=false/true
molecular=linear/circle
number3errors=0
minlen=100
maxlen=5000 
FRpairs=false/true
SequenceExtract=false/true
CTconversion=false/true
ShowPrimerAlignment=true/false
ShowOnlyAmplicons=true/false
ShowPCRProducts=true/false
ShowPrimerAlignment=true/false
ShowPrimerAlignmentPCRproduct=true/false 
ShowPCRproductCalculation=true/false
```

 
**minlen=/maxlen=**
> The box of “Minimal and Maximal PCR Product length (bp)” – has the default value of 5000 bp, allowing the user to define the maximal size of the expected PCR product. Any amplicons larger than a defined value will be filtered out. 

**ShowPCRProducts=true/false**
> “PCR product prediction” has the default value checked; to search for primer binding sites without further analysis of potential PCR fragments, this option should be disabled.

**molecular=linear/circle**
>“Circular sequence” – analysis of circular molecules (plasmid, mitochondria or plastids DNA, etc.); in this case, the primers can produce one or two amplicons.

**CTconversion=false/true**
> “C >> T bisulphite conversion” – bisulphite modified genome sequence, design of specific PCR primers for in silico bisulphite conversion for both strands - only cytosines not followed by guanidine (CpG methylation) will be replaced by thymines:
```
5’aaCGaagtCCCCa-3'        5’aaCGaagtTTTTa-3'
  |||||||||||||     ->      ||||||:|::::| 
3’ttGCttCaggggt-5'        3’ttGCttTaggggt
```
**FRpairs=false/true**
> “Restrict analysis for F/R primer pairs” – this option is used to analyse the primers list, where the common name unites each primer pair. Moreover, the program is not limited to one unique pair per primer: for one "Forward" primer, there can be several "reverse" primers, the same as the “reverse” primer. The search for potential amplicons will be carried out only for these primer pairs, while other primers from the list will be ignored. It is possible to group pairs in the form of a table, with each line containing a separate pair:
```
1020	aggcctgtgatgctgatgat	cccaacaccaaagaggaaag
1021	gctctgacctattgtctgtctgtct	cagtctccacagcagcagag
1022	gtcctgctgacacacaccact	cgcaggagtagaagaaaga
1023	gatagggaggtgggggtct	tgattggacaggcagcacag
1024	catgtgagagggagggcta	gggactgcatgctggtg
1025	cacaacgagagctggggaga	ggataattgctgcaagagagaaa
1026	ccttggagggggatgtagag	tactacggccgcgagg
1027	gaaggatctttacccctctctcc	tcccagggtgcggagg
```

**ShowPrimerAlignment=true/false**
> “Show all matching for primers alignment” - checked by default, the software shows the result, including all matching of stable binding primer to the target. Not in all cases combinations of primers can produce the PCR products in the current assay conditions. Still, the user can examine the stability of primer binding sites, orientation and coordinates in the target.

**ShowPrimerAlignmentPCRproduct=true/false** 
> “Show alignment only for matching primers for PCR product” - all primer binding sites were represented in the previous option. In this case, the primer and target alignment analysis will be shown only for matching primers.

**ShowOnlyAmplicons=true/false**
“Show only amplicon lengths” – checking this option allows the user to collect only amplicon lengths without primer and target alignment analysis. This option is recommended for in silico PCR of the whole genome, including all chromosome analysis with highly abandoned repeated sequences (in silico PCR for techniques based on repeats: iPBS, IRAP, ISSR or RAPD).

**linkedsearch=false/true**
> “Linked (Associated) search” - programmable searching when binding sites for primers or probes are searched within a determined distance.  In linked searching, the search criteria are based only on the distance between forward sites. Linked searching can perform a broader variety of tasks, ranging from conventional sequence matching to in silico PCR and general tasks involving DNA sequence analysis with approximate matching.
> To better understand how linked searching works in practice, we use the example of in silico PCR with two degenerate Copia-type RT primers. The primer sequences were converted into a single line, where the forward primer sequence (5′-CARATGGAYGTNAARAC) was followed by the expected distance between the primer-binding sites (200–300 nt) and the complementary sequence (TAYGTNGAYGAYATG) of the reverse primer (5′-CATRTCRTCNACRTA):  
```
>RT+(QMDVK)_RT-(YVDDML)
CARATGGAYGTNAARAC[200-300]TAYGTNGAYGAYATG

or without specifying the interval between sequences:

>RT+(QMDVK)_RT-(YVDDML)
CARATGGAYGTNAARAC[300]TAYGTNGAYGAYATG

or if it is specified that the second sequence should be made complementary:

>RT+(QMDVK)_RT-(YVDDML)
CARATGGAYGTNAARAC[300]@CATRTCRTCNACRTA
```

**type=primer/probe**
> “Probe search” – helps the user execute searching of binding sites not only for primers but also for probes (TaqMan, Molecular Beacon, microRNA (miRNA), microarrays, etc.). This option is recommended in the cases when primer binding sites were not found or for searching for binding sites of probes for which the complementarity is expected only for part of the sequence, for example, in “Molecular Beacon” (both termini have not complementary regions to the target).

 


## Sequence Entry:

Sequence data files are prepared using a text editor and saved in ASCII as text/plain format (.txt) or in .fasta or without file extensions (a file extension is not obligatory). The program takes a single sequence or accepts multiple DNA/RNA sequences in FASTA format. The template length is not limited.

## FASTA format description:
A sequence in FASTA format consists of:
One line starts with a ">" sign and a sequence identification code. A textual description of the sequence optionally follows it. Since it is not part of the official format description, the software can ignore it when it is present.
One or more lines containing the sequence itself. A file in FASTA format may comprise more than one sequence.
The input DNA/RNA sequence can contain the degenerate nucleotides accepted as IUPAC code, an extended vocabulary of 13 letters, which allows the description of ambiguous DNA code. Each letter represents a combination of one or several nucleotides: M (A/C), R (A/G), W (A/T), S (G/C), Y (C/T), K (G/T), V (A/G/C), H (A/C/T), D (A/G/T), B (C/G/T), N (A/G/C/T), U (Uracil), I (Inosine). LNA: dA=E, dC=F, dG=J, dT=L.


## The output is saved in tab-delimited, plain text files. 


## In silico PCR application, primers file examples:

For the beginning, the list of query primer(s) should be prepared in the FASTA format or as a table (with TAB or whitespace separators) (table from Microsoft Office documents or Excel worksheet) or only bare sequence without space between letters. The software reads only the first two columns with names and sequences or the first three/four columns for primers and probes:

```
ITS1	TCCGTAGGTGAACCTGCGG
ITS2	GCTGCGTTCTTCATCGATGC
ITS3	GCATCGATGAAGAACGCAGC
ITS4	TCCTCCGCTTATTGATATGC
ITS5	GGAAGTAAAAGTCGTAACAAGG
KAN2-FP	ACCTACAACAAAGCTCTCATCAACC
KAN2-RP	GCAATGTAACATCAGAGATTTTGAG
L_SP6	TCAAGCTATGCATCCAACGCG
L_T7	TAGGGCGAATTGGGCCCGACG
```
The first column indicates the primer's name, while the second column contains the primer sequence. This is the most convenient format for storing and using primers in such studies. Within primers, sequence the spaces and no DNA letters are allowed. Primer names can contain any characters, including only space. Furthermore, the names of primers can be identical.
FASTA format has a description line starting with the “>” sign followed by a plain DNA sequence. This format is widespread for the storage and processing of DNA sequences:
```
>ITS1 
TCCGTAGGTGAACCTGCGG
>ITS2 
GCTGCGTTCTTCATCGATGC
>ITS3 
GCATCGATGAAGAACGCAGC
>ITS4 
TCCTCCGCTTATTGATATGC
>ITS5 
GGAAGTAAAAGTCGTAACAAGG
>KAN2-FP 
ACCTACAACAAAGCTCTCATCAACC
>KAN2-RP 
GCAATGTAACATCAGAGATTTTGAG
>L_SP6 
TCAAGCTATGCATCCAACGCG
>L_T7 
TAGGGCGAATTGGGCCCGACG
```


