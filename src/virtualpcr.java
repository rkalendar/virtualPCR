import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.InvalidPathException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class virtualpcr {

    public static void main(String[] args) {
        if (args.length == 0 || isHelpRequested(args)) {
            printUsage();
            return;
        }
        if (args.length > 0) {
            String tagfile = "";
            String outfile = "";
            String primersfile = "";
            String infile = args[0];

            int minlen = 30;
            int maxlen = 3000;
            int Err3end = 1;

            boolean isprobe = false;
            boolean ispattern = false;
            boolean iscircle = false;
            boolean seqextract = false;
            boolean FRpairs = false;                 //SetLookR_Fprimes
            boolean pcr_predict = true;              //SetShowPCRProducts
            boolean alignment = true;                //SetShowPrimerAlignment
            boolean ShowOnlyAmplicons = false;       //SetShowOnlyAmplicons
            boolean PCRmatch_alignment = true;       //SetShowPrimerAlignmentPCRproduct
            boolean CTconversion = false;            //CTconversion=false/true
            boolean PrimerStatistic = true;          //Statistic report for each primer

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            try (BufferedReader br = new BufferedReader(new FileReader(infile))) {
                String line;

                while ((line = br.readLine()) != null) {
                    String cline = line;
                    line = line.toLowerCase();

                    if (line.contains("output_path=")) {
                        String value = line.substring(line.indexOf('=') + 1).trim(); // safer than hardcoded 13
                        if (!value.isEmpty()) {
                            outfile = value;
                        }
                    }

                    if (line.contains("targets_path=")) {
                        tagfile = cline.substring(line.indexOf('=') + 1).trim(); // tolerate leading whitespace
                        System.out.println(tagfile);
                    }
                    if (line.contains("primers_path=")) {
                        primersfile = cline.substring(line.indexOf('=') + 1).trim();
                        System.out.println(primersfile);
                    }
                    if (line.contains("primerstatistic=false")) {
                        PrimerStatistic = false;
                        System.out.println(line);
                    }
                    if (line.contains("type=probe")) {
                        isprobe = true;
                        System.out.println(line);
                    }
                    if (line.contains("molecular=circle")) {
                        iscircle = true;
                        System.out.println(line);
                    }
                    if (line.contains("linkedsearch=true")) {
                        ispattern = true;
                        System.out.println(line);
                    }
                    if (line.contains("showprimeralignmentpcrproduct=true")) {
                        PCRmatch_alignment = true;
                        System.out.println(line);
                    }
                    if (line.contains("sequenceextract=true")) {
                        seqextract = true;
                        System.out.println(line);
                    }
                    if (line.contains("frpairs=true")) {
                        FRpairs = true;
                        System.out.println(line);
                    }
                    if (line.contains("showprimeralignment=false")) {
                        alignment = false;
                        System.out.println(line);
                    }
                    if (line.contains("showpcrproducts=false")) {
                        pcr_predict = false;
                        System.out.println(line);
                    }

                    if (line.contains("showonlyamplicons=true")) {
                        ShowOnlyAmplicons = true;
                        System.out.println(line);
                    }
                    if (line.contains("ctconversion=true")) {
                        CTconversion = true;
                        System.out.println(line);
                    }
                    if (line.contains("number3errors=")) {
                        int h = StrToInt(line.substring(line.indexOf('=') + 1));
                        System.out.println(line);
                        if (h < 0) {
                            h = 0;
                        }
                        if (h > 10) {
                            h = 10;
                        }
                        Err3end = h;
                    }
                    if (line.contains("minlen=")) {
                        int h = StrToInt(line.substring(line.indexOf('=') + 1));
                        System.out.println(line);
                        if (h < 20) {
                            h = 20;
                        }
                        if (h > 50000) {
                            h = 50000;
                        }
                        minlen = h;
                    }
                    if (line.contains("maxlen=")) {
                        int h = StrToInt(line.substring(line.indexOf('=') + 1));
                        System.out.println(line);
                        if (h < minlen) {
                            h = minlen;
                        }
                        if (h > 50000) {
                            h = 50000;
                        }
                        maxlen = h;
                    }
                }
            } catch (IOException e) {
                System.out.println("Failed to read config file: " + infile + " - " + e.getMessage());
                return;
            }

            System.out.println("\n");

            PrimerLoader.Primers pr = PrimerLoader.loadPrimers(primersfile);
            String[] PrimersList = pr.list;
            String[] PrimersName = pr.names;
            String[] PrimersOri = pr.originals;
            int n = pr.count;

            if (n > 0) {
                File folder = new File(tagfile);
                if (folder.exists() && (folder.isDirectory() || folder.isFile())) {

                    if (outfile.isBlank()) {
                        Path p = folder.toPath().toAbsolutePath().normalize();
                        Path outDir = Files.isDirectory(p) ? p : (p.getParent() != null ? p.getParent() : p);
                        outfile = outDir.toString() + File.separator + "report.out";
                    } else {
                        Path outPath;
                        try {
                            OutputPath.Target tgt = OutputPath.parse(outfile);
                            OutputPath.ensureParentExists(tgt);
                            outPath = (tgt.type() == OutputPath.Type.DIRECTORY)
                                    ? tgt.path().resolve("report.out")
                                    : tgt.path();
                        } catch (InvalidPathException | IOException e) {
                            outPath = folder.toPath().resolve("report.out");
                        }
                        outfile = outPath.toAbsolutePath().toString();
                    }
                    System.out.println("\noutput_path=" + outfile + "\n");

                    if (folder.isDirectory()) {
                        File[] files = folder.listFiles();
                        if (files == null) { // unreadable directory / I/O error / race after isDirectory()
                            System.out.println("\nFailed to read directory: " + folder);
                            return;
                        }
                        int k = -1;
                        String[] filelist = new String[files.length];
                        for (File file : files) {
                            if (file.isFile()) {
                                filelist[++k] = file.getAbsolutePath();
                            }
                        }

                        if (PrimerStatistic) {
                            InSilicoPCR.beginGlobalStat();
                        }
                        // Single UTF-8 writer owns the combined report: one open, one copy of each
                        // file's section, then the global summary. Avoids the previous double-write.
                        try {
                            Path outPath = Paths.get(outfile);
                            Path parent = outPath.getParent();
                            if (parent != null) {
                                Files.createDirectories(parent);
                            }
                            int processed = 0;
                            try (BufferedWriter w = Files.newBufferedWriter(outPath, StandardCharsets.UTF_8)) {
                                for (String nfile : filelist) {
                                    if (nfile == null) {
                                        continue; // listFiles() may leave trailing nulls for sub-directories
                                    }
                                    try {
                                        StringBuilder sr = Run(nfile, PrimersList, PrimersName, PrimersOri, isprobe, ispattern, iscircle, seqextract, Err3end, minlen, maxlen, FRpairs, pcr_predict, alignment, ShowOnlyAmplicons, PCRmatch_alignment, CTconversion, PrimerStatistic);
                                        if (sr == null) {
                                            System.out.println("Skipped (no sequences): " + nfile);
                                            continue;
                                        }
                                        w.write(sr.toString());
                                        w.write("\n\n");
                                        processed++;
                                        System.out.println("Processed: " + nfile);
                                    } catch (RuntimeException ex) {
                                        // Isolate one malformed file so it cannot abort the whole batch.
                                        System.out.println("Skipped (error): " + nfile + " - " + ex);
                                    }
                                }
                                if (PrimerStatistic) {
                                    String summary = InSilicoPCR.getGlobalStatReport();
                                    if (!summary.isEmpty()) {
                                        w.write(summary);
                                        System.out.print(summary);
                                    }
                                }
                            }
                            if (processed == 0) {
                                System.out.println("Warning: no target files produced output in " + folder);
                            }
                            System.out.println("Saved report: " + outfile);
                            if (PrimerStatistic) {
                                writeStatsTsv(outfile);
                            }
                        } catch (IOException e) {
                            System.out.println("Failed to write report: " + e.getMessage());
                        }

                    } else {
                        if (PrimerStatistic) {
                            InSilicoPCR.beginGlobalStat();
                        }
                        StringBuilder sr = Run(tagfile, PrimersList, PrimersName, PrimersOri, isprobe, ispattern, iscircle, seqextract, Err3end, minlen, maxlen, FRpairs, pcr_predict, alignment, ShowOnlyAmplicons, PCRmatch_alignment, CTconversion, PrimerStatistic);
                        if (sr == null) {
                            System.out.println("No target sequences were processed in: " + tagfile);
                        } else {
                            if (PrimerStatistic) {
                                String summary = InSilicoPCR.getGlobalStatReport();
                                if (!summary.isEmpty()) {
                                    sr.append(summary);
                                }
                            }
                            System.out.println(sr);
                            try {
                                writeFile(outfile, sr.toString());
                                System.out.println("Saved report: " + outfile);
                            } catch (IOException e) {
                                System.out.println("Failed to write report: " + e.getMessage());
                            }
                            if (PrimerStatistic) {
                                writeStatsTsv(outfile);
                            }
                        }
                    }
                } else {
                    System.out.println("\nFailed to open file: " + folder);
                }
            } else {
                System.out.println("\nFailed to open primer's file: " + primersfile);
            }
        }
    }

    private static boolean isHelpRequested(String[] args) {
        for (String arg : args) {
            if (arg == null) {
                continue;
            }
            switch (arg.trim().toLowerCase()) {
                case "-help", "--help", "-h", "--h", "/help", "/h", "/?", "-?", "?", "help" -> {
                    return true;
                }
                default -> {
                }
            }
        }
        return false;
    }

    private static void printUsage() {
        System.out.print("""
            virtualPCR (2024-2026) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)
            https://github.com/rkalendar/virtualPCR
            In silico PCR for simple and complex tasks.

            USAGE:
              java -jar virtualPCR.jar <config.file>
              java -jar virtualPCR.jar -help

              For large genomes, allocate more heap memory:
              java -Xms4g -Xmx16g -jar virtualPCR.jar <config.file>

            HELP:
              -help, --help, -h, /help, /?, -?, ?   Show this help and exit.
              Running with no arguments also shows this help.

            CONFIGURATION FILE:
              All parameters are given in a plain-text file, one key=value per line.
              Option keys are case-insensitive; boolean options accept true or false.

              PATHS
                targets_path=PATH   Target sequence file, or a directory of files (FASTA or text).
                primers_path=PATH   Primer/probe file (FASTA, tab- or space-delimited).
                output_path=PATH    Output file or directory; empty = write next to the target(s).

              SEARCH MODE
                type=primer|probe              Search mode (default: primer).
                molecular=linear|circle        Template topology (default: linear).
                linkedsearch=true|false        Linked/associated search (default: false).
                frpairs=true|false             Restrict to defined F/R primer pairs (default: false).
                ctconversion=true|false        Simulate bisulfite C->T conversion (default: false).

              FILTERS
                minlen=N                       Min amplicon length, bp; 20-50000 (default: 30).
                maxlen=N                       Max amplicon length, bp; minlen-50000 (default: 3000).
                number3errors=N                Mismatches allowed near the 3' end; 0-10 (default: 1).

              OUTPUT CONTROL
                primerstatistic=true|false                Per-primer summary statistics (default: true).
                sequenceextract=true|false                Extract amplicon sequences (default: false).
                showprimeralignment=true|false            Show primer-target alignments (default: true).
                showpcrproducts=true|false                Report predicted PCR products (default: true).
                showonlyamplicons=true|false              Report amplicon lengths only (default: false).
                showprimeralignmentpcrproduct=true|false  Only product-forming alignments (default: true).

            EXAMPLE config.file:
                targets_path=test/ch02.fasta
                primers_path=test/its.txt
                type=primer
                molecular=linear
                minlen=100
                maxlen=2000
                number3errors=1

            INPUT FORMATS:
              Target : plain text, or single-/multi-entry FASTA; length limited only by memory.
              Primers: FASTA, or tab/space-delimited (col 1 = ID, remaining cols = sequence(s)).
                       Sequences < 12 nt are skipped. IUPAC, LNA (E,F,J,L) and inosine (I) supported.

            OUTPUT:
              Tab-delimited plain-text report: binding-site coordinates, alignments, predicted
              product sizes, extracted amplicons, and per-primer statistics (per options enabled).

            Documentation:  https://github.com/rkalendar/virtualPCR
            Online version: https://primerdigital.com/tools/epcr.html
            """);
    }

    // Computes the report for one target file and returns it. The caller owns all file output.
    private static StringBuilder Run(String tagfile, String[] PrimersList, String[] PrimersName, String[] PrimersOri, boolean isprobe, boolean ispattern, boolean iscircle, boolean seqextract, int Err3end, int minlen, int maxlen, boolean FRpairs, boolean pcr_predict, boolean alignment, boolean ShowOnlyAmplicons, boolean PCRmatch_alignment, boolean CTconversion, boolean PrimerStatistic) {
        StringBuilder sr = new StringBuilder(1000);
        try {
            System.out.println("Running...");
            System.out.println("\nTarget file name: " + tagfile);

            InSilicoPCR s2 = new InSilicoPCR();
            ReadingSequencesFiles rf = new ReadingSequencesFiles(Paths.get(tagfile));
            if (rf.getNseq() == 0) {
                System.out.println("There is no sequence(s).");
                System.out.println("File format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
                System.out.println(">seq1\nactacatactacatcactctctctccgcacag\n");
                return null; // no sequences -> signal "not processed" to the caller
            }
            System.out.println("Target sequence length = " + rf.getLength() + " nt");
            s2.SetSequences(rf.getSequences(), rf.getNames());
            s2.SetCurrentFileName(tagfile);
            s2.SetPrimerStat(PrimerStatistic);
            s2.SetLookR_Fprimes(FRpairs);
            s2.SetShowPCRProducts(pcr_predict);
            s2.SetShowPrimerAlignment(alignment);
            s2.SetShowOnlyAmplicons(ShowOnlyAmplicons);
            s2.SetShowPrimerAlignmentPCRproduct(PCRmatch_alignment);
            s2.SetProductMaxLength(maxlen);
            s2.SetProductMinLength(minlen);
            s2.SetCTbisulfate(CTconversion);
            s2.SetPrimers(PrimersList, PrimersName, PrimersOri, isprobe, ispattern, iscircle, seqextract, Err3end);
            long startTime = System.nanoTime();
            s2.Run();

            long duration = (System.nanoTime() - startTime) / 1000000;
            StringBuilder st = new StringBuilder();
            if (duration > 999) {
                duration = duration / 1000;
                st.append("Time taken: ").append(duration).append(" seconds\n\n");
            } else {
                st.append("Time taken: ").append(duration).append(" milliseconds\n\n");
            }

            sr.append(s2.getResult()).append(st);

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
            return null; // load failure -> signal "not processed" to the caller
        }
        return sr;
    }

    // Write text to a file as UTF-8, creating parent directories and truncating any existing file.
    private static void writeFile(String outfile, CharSequence text) throws IOException {
        Path path = Paths.get(outfile);
        Path parent = path.getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        try (BufferedWriter w = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            w.append(text);
        }
    }

    // Derive the stats-file path from the report path: report.out -> report.stats.tsv
    private static String deriveStatsPath(String outfile) {
        int dot = outfile.lastIndexOf('.');
        int sep = Math.max(outfile.lastIndexOf('/'), outfile.lastIndexOf('\\'));
        return (dot > sep) ? outfile.substring(0, dot) + ".stats.tsv" : outfile + ".stats.tsv";
    }

    // Write the global per-primer statistics as a stand-alone TSV next to the report (if any).
    private static void writeStatsTsv(String outfile) {
        String tsv = InSilicoPCR.getGlobalStatTsv();
        if (tsv.isEmpty()) {
            return;
        }
        String statsFile = deriveStatsPath(outfile);
        try {
            writeFile(statsFile, tsv);
            System.out.println("Saved statistics: " + statsFile);
        } catch (IOException e) {
            System.out.println("Failed to write statistics TSV.");
        }
    }

    public static int StrToInt(String str) {
        StringBuilder r = new StringBuilder();
        int z = 0;
        r.append(0);
        for (int i = 0; i < str.length(); i++) {
            char chr = str.charAt(i);
            if (chr > 47 && chr < 58) {
                r.append(chr);
                z++;
                if (z > 10) {
                    break;
                }
            }
            if (chr == '.' || chr == ',') {
                break;
            }
        }
        long v = Long.parseLong(r.toString()); // up to 11 digits fits in long; avoids int-overflow NumberFormatException
        return (v > Integer.MAX_VALUE) ? Integer.MAX_VALUE : (int) v;
    }
}
