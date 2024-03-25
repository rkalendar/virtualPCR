## virtualPCR v1.0
## "in silico PCR for simple and complex tasks"
by Ruslan Kalendar 

email: ruslan.kalendar@helsinki.fi

[Web](https://primerdigital.com/tools/)

## Availability and requirements:

Operating system(s): Platform independent

Programming language: Java 20 or higher

[Java Downloads](https://www.oracle.com/java/technologies/downloads/)


How do I set or change [the Java path system variable](https://www.java.com/en/download/help/path.html)


To run the project from the command line, go to the target folder and type the following; an individual file or a file folder can be specified:

```java -jar \virtualPCR\dist\virtualPCR.jar \virtualPCR\task\config.file ```

### Basic usage: 
To enter parameters and specify the location of the target files and primer file, you must specify this via a file on the command line. An example of such a file here (file name or extension does not matter):

> **config.file**
```
targets_path=C:\virtualPCR\test\1.txt
primers_path=C:\virtualPCR\task\primers.txt
type=primer/probe
linkedsearch=false/true
molecular=linear/circle
number3errors=1
minlen=100
maxlen=5000 
FRpairs=false
ShowPrimerAlignment=true/false
ShowOnlyAmplicons=true/false
ShowPCRProducts=true/false
ShowPrimerAlignment=true/false
ShowPrimerAlignmentPCRproduct=true/false 
ShowPCRproductCalculation=true/false
```

> -maxlen=/maxlen=
> The box of “Minimal and Maximal PCR Product length (bp)” – has the default value of 5000 bp, allowing the user to define the maximal size of the expected PCR product. Any amplicons larger than a defined value will be filtered out. 

>-ShowPCRProducts=true -“PCR product prediction” – has the default value of checked; to search for primer binding sites without further analysis of potential PCR fragments, this option should be disabled.
-molecular=linear/circle - “Circular sequence” – analysis of circular molecules (plasmid, mitochondria or plastids DNA, etc.); in this case, the primers can produce one or two amplicons.

> -CTconversion=false/true
> “C >> T bisulphite conversion” ¬– bisulphite modified genome sequence, design of specific PCR primers for in silico bisulphite conversion for both strands - only cytosines not followed by guanidine (CpG methylation) will be replaced by thymines.

> FRpairs=false
> “Restrict analysis for F/R primer pairs” – this option is used to analyse the primers list, where the common name unites each primer pair. A similar analysis can be carried out for primers with the same names or names that differ in the last letter – F/R. The program will recognize paired primers (Forward as F and Reverse as R). For this type of analysis, primers from the same pair(s) must have identical names but finish using “R or F” (e.g. “seq1R” and “seq1F” form a pair). The name length and structure (including "F" and "R" inside names) are not important. Moreover, the program is not limited to one unique pair per primer: for one "Forward" primer, there can be several "reverse" primers, the same as the “reverse” primer. The search for potential amplicons will be carried out only for these primer pairs, while other primers from the list will be ignored.



## Sequence Entry:

Sequence data files are prepared using a text editor and saved in ASCII as text/plain format (.txt) or in .fasta or without file extensions (a file extension is not obligatory). The program takes a single sequence or accepts multiple DNA sequences in FASTA format. The template length is not limited.

## FASTA format description:
A sequence in FASTA format consists of:
One line starts with a ">" sign and a sequence identification code. A textual description of the sequence optionally follows it. Since it is not part of the official format description, the software can ignore it when it is present.
One or more lines containing the sequence itself. A file in FASTA format may comprise more than one sequence.
The input DNA sequence can contain the degenerate nucleotides accepted as IUPAC code, an extended vocabulary of 11 letters, which allows the description of ambiguous DNA code. Each letter represents a combination of one or several nucleotides: M (A/C), R (A/G), W (A/T), S (G/C), Y (C/T), K (G/T), V (A/G/C), H (A/C/T), D (A/G/T), B (C/G/T), N (A/G/C/T), and also U (T) and I (Inosine).


## The output is saved in tab-delimited, plain text files. 
