## virtualPCR
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
To enter parameters and specify the location of the target files and primer file, you must specify this via file on the command line. An example of such a file, here (file name or extension, does not matter):

> **config.file**
```
targets_path=C:\virtualPCR\test\1.txt
primers_path=C:\virtualPCR\task\primers.txt
type=primer/probe
linkedsearch=false/true
molecular=linear/circle
number3errors=1
minlen=100
maxlen=5500
FRpairs=false
ShowPrimerAlignment=true/false
ShowOnlyAmplicons=true/false
ShowPCRProducts=true/false
ShowPrimerAlignment=true/false
ShowPrimerAlignmentPCRproduct=true/false 
ShowPCRproductCalculation=true/false
```

## Sequence Entry:

Sequence data files are prepared using a text editor and saved in ASCII as text/plain format (.txt) or in .fasta or without file extensions (a file extension is not obligatory). The program takes a single sequence or accepts multiple DNA sequences in FASTA format. The template length is not limited.

## FASTA format description:
A sequence in FASTA format consists of:
One line starts with a ">" sign and a sequence identification code. A textual description of the sequence optionally follows it. Since it is not part of the official format description, software can ignore it when it is present.
One or more lines containing the sequence itself. A file in FASTA format may comprise more than one sequence.



## The output is saved in file: tab-delimited, plain text files. 
