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
        {
            String tagfile = "";
            String outfile = "";
            String primersfile = "";
            String infile = args[0];
            for (int i = 1; i < args.length; i++) { // only the config file is read; say so rather than ignoring it
                System.out.println("Warning: argument ignored (only the config file is used): " + args[i]);
            }

            int minlen = 30;
            int maxlen = 3000;
            int Err3end = 1;
            int flanks = 0;

            boolean isprobe = false;
            boolean ispattern = false;
            boolean iscircle = false;
            boolean seqextract = false;
            boolean FRpairs = false;                 //SetLookR_Fprimes
            boolean pcr_predict = true;              //SetShowPCRProducts
            boolean alignment = true;                //SetShowPrimerAlignment
            boolean ShowOnlyAmplicons = false;       //SetShowOnlyAmplicons
            boolean PCRmatch_alignment = true;       //SetShowPrimerAlignmentPCRproduct
            boolean PrimerStatistic = true;          //Statistic report for each primer

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            try (BufferedReader br = new BufferedReader(new FileReader(infile))) {
                String line;

                while ((line = br.readLine()) != null) {
                    String echo = line.toLowerCase();          // console log kept as in earlier versions
                    String raw = stripComment(line).trim();
                    if (raw.isEmpty()) {
                        continue;
                    }
                    int eq = raw.indexOf('=');
                    if (eq < 0) {
                        continue;                              // a note, a command line, a bare path - not an option
                    }
                    String key = raw.substring(0, eq).trim().toLowerCase();
                    String value = raw.substring(eq + 1).trim();   // original case: paths are case-sensitive
                    String lval = value.toLowerCase();
                    if (!key.matches("[a-z0-9_]+")) {
                        continue;                              // whatever this line is, it is not an option either
                    }
                    if (!KNOWN.contains(key)) {
                        // Previously such a line fell through the whole if-chain without a word, so a
                        // typo like min_len= or primer_path= silently ran the whole job on defaults.
                        System.out.println("Warning: unknown option ignored: " + key + "=" + value);
                        continue;
                    }

                    switch (key) {
                        case "output_path":
                            if (!value.isEmpty()) {
                                outfile = value;
                            }
                            break;
                        case "targets_path":
                            tagfile = value;
                            System.out.println(tagfile);
                            break;
                        case "primers_path":
                            primersfile = value;
                            System.out.println(primersfile);
                            break;
                        case "primerstatistic":
                            PrimerStatistic = !lval.startsWith("false");
                            System.out.println(echo);
                            break;
                        case "type":
                            isprobe = lval.startsWith("probe");
                            System.out.println(echo);
                            break;
                        case "molecular":
                            iscircle = lval.startsWith("circle");
                            System.out.println(echo);
                            break;
                        case "linkedsearch":
                            ispattern = lval.startsWith("true");
                            System.out.println(echo);
                            break;
                        case "showprimeralignmentpcrproduct":
                            // Defaults to true, so anything but an explicit "false" leaves it on.
                            PCRmatch_alignment = !lval.startsWith("false");
                            System.out.println(echo);
                            break;
                        case "sequenceextract":
                            seqextract = lval.startsWith("true");
                            System.out.println(echo);
                            break;
                        case "flanks":
                        case "flanking":
                            // flanks=N: N extra bases of the template on each side of every extracted amplicon
                            System.out.println(echo);
                            flanks = readInt(key, value, 0, 50000);
                            break;
                        case "frpairs":
                            FRpairs = lval.startsWith("true");
                            System.out.println(echo);
                            break;
                        case "showprimeralignment":
                            alignment = !lval.startsWith("false");
                            System.out.println(echo);
                            break;
                        case "showpcrproducts":
                            pcr_predict = !lval.startsWith("false");
                            System.out.println(echo);
                            break;
                        case "showonlyamplicons":
                            ShowOnlyAmplicons = lval.startsWith("true");
                            System.out.println(echo);
                            break;
                        case "number3errors":
                            System.out.println(echo);
                            Err3end = readInt(key, value, 0, 10);
                            break;
                        case "minlen":
                            System.out.println(echo);
                            minlen = readInt(key, value, 20, 50000);
                            break;
                        case "maxlen":
                            // Not compared with minlen here: these are separate lines of one file, so a
                            // "maxlen=" preceding "minlen=" would be checked against the default. The two
                            // are reconciled once, after the whole file is read.
                            System.out.println(echo);
                            maxlen = readInt(key, value, 0, 50000);
                            break;
                        default:
                            break;
                    }
                }
            } catch (IOException e) {
                System.out.println("Failed to read config file: " + infile + " - " + e.getMessage());
                return;
            }

            if (maxlen < minlen) { // reconciled here, so the order of the two lines in the file does not matter
                System.out.println("maxlen=" + maxlen + " < minlen=" + minlen + " -> maxlen=" + minlen);
                maxlen = minlen;
            }

            if (flanks > 0 && !seqextract) { // flanking regions are meaningless without the sequence they flank
                seqextract = true;
                System.out.println("flanks=" + flanks + " -> sequenceextract=true (enabled automatically)");
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
                            // resolve() on a plain file would build ".../ch02.fasta/report.out";
                            // fall back to the directory holding the target, as the blank branch does.
                            System.out.println("Unusable output_path=" + outfile + " (" + e.getMessage()
                                    + "); writing next to the target(s) instead.");
                            Path p = folder.toPath().toAbsolutePath().normalize();
                            Path outDir = Files.isDirectory(p) ? p : (p.getParent() != null ? p.getParent() : p);
                            outPath = outDir.resolve("report.out");
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
                                        StringBuilder sr = Run(nfile, PrimersList, PrimersName, PrimersOri, isprobe, ispattern, iscircle, seqextract, flanks, Err3end, minlen, maxlen, FRpairs, pcr_predict, alignment, ShowOnlyAmplicons, PCRmatch_alignment, PrimerStatistic);
                                        if (sr == null) {
                                            System.out.println("Skipped (no sequences): " + nfile);
                                            continue;
                                        }
                                        writeChunked(w, sr);
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
                        StringBuilder sr = Run(tagfile, PrimersList, PrimersName, PrimersOri, isprobe, ispattern, iscircle, seqextract, flanks, Err3end, minlen, maxlen, FRpairs, pcr_predict, alignment, ShowOnlyAmplicons, PCRmatch_alignment, PrimerStatistic);
                        if (sr == null) {
                            System.out.println("No target sequences were processed in: " + tagfile);
                        } else {
                            if (PrimerStatistic) {
                                String summary = InSilicoPCR.getGlobalStatReport();
                                if (!summary.isEmpty()) {
                                    sr.append(summary);
                                }
                            }
                            // Printing the whole report meant another full copy plus minutes of
                            // scrolling on a genome-sized run. The head keeps short runs unchanged
                            // and still shows something for a piped stdout.
                            if (sr.length() > 8192) {
                                System.out.println(sr.substring(0, 8192));
                                System.out.println("...[truncated, see " + outfile + "]\n");
                            } else {
                                System.out.println(sr);
                            }
                            try {
                                writeFile(outfile, sr);
                                System.out.println("Saved report: " + outfile);
                            } catch (IOException e) {
                                System.out.println("Failed to write report: " + e.getMessage());
                            }
                            if (PrimerStatistic) {
                                writeStatsTsv(outfile);
                            }
                        }
                    }
                } else if (tagfile.isBlank()) {
                    System.out.println("\nNo targets_path= in the config file: " + infile);
                } else {
                    System.out.println("\nFailed to open file: " + folder);
                }
            } else if (primersfile.isBlank()) {
                System.out.println("\nNo primers_path= in the config file: " + infile);
            } else if (!new File(primersfile).isFile()) {
                // PrimerLoader already printed this exact sentence; repeat it only when the file
                // really is missing, so the wording keeps matching the cause. Greps rely on it.
                System.out.println("\nFailed to open primer's file: " + primersfile);
            } else {
                System.out.println("\nNo usable primers in " + primersfile
                        + " (need ID + sequence per line, or FASTA; processed sequences of 11 characters or fewer are skipped).");
            }
        }
    }

    // Every option the parser acts on. A line whose key is not here is echoed as a warning
    // instead of falling through unnoticed. Both spellings of the flank option are accepted.
    private static final java.util.Set<String> KNOWN = java.util.Set.of(
            "output_path", "targets_path", "primers_path", "primerstatistic", "type", "molecular",
            "linkedsearch", "showprimeralignmentpcrproduct", "sequenceextract", "flanks", "flanking",
            "frpairs", "showprimeralignment", "showpcrproducts", "showonlyamplicons",
            "number3errors", "minlen", "maxlen");

    // A '#' is legal inside a path, so it only opens a comment where it cannot belong to a value:
    // before the '=' of a key=value pair, or set off by whitespace.
    private static String stripComment(String s) {
        int eq = s.indexOf('=');
        for (int i = 0; i < s.length(); i++) {
            if (s.charAt(i) == '#'
                    && ((eq < 0 || i < eq) || (i > 0 && Character.isWhitespace(s.charAt(i - 1))))) {
                return s.substring(0, i);
            }
        }
        return s;
    }

    // StrToInt cannot report failure - it just concatenates the digits it meets - so anything the
    // user wrote that is not a plain number, and any value the limits had to move, is reported here.
    // Deliberately not a hard filter: "minlen=1000bp" and "flanks=100 nt" keep working as before.
    private static int readInt(String key, String value, int lo, int hi) {
        int h = StrToInt(value);
        if (!value.matches("[0-9]+")) {
            System.out.println("Note: " + key + "=" + value + " read as " + h);
        }
        int c = (h < lo) ? lo : ((h > hi) ? hi : h);
        if (c != h) {
            System.out.println("Note: " + key + "=" + h + " outside " + lo + ".." + hi + " -> " + c);
        }
        return c;
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

              FILTERS
                minlen=N                       Min amplicon length, bp; 20-50000 (default: 30).
                maxlen=N                       Max amplicon length, bp; minlen-50000 (default: 3000).
                number3errors=N                Mismatches allowed near the 3' end; 0-10 (default: 1).

              OUTPUT CONTROL
                primerstatistic=true|false                Per-primer summary statistics (default: true).
                sequenceextract=true|false                Extract amplicon sequences (default: false).
                flanks=N                                  Flanking bases kept on each side of the amplicon;
                                                          0-50000 (default: 0); implies sequenceextract=true.
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
                sequenceextract=true
                flanks=100

            INPUT FORMATS:
              Target : single-/multi-entry FASTA (each entry must start with a '>' header line);
                       length limited only by memory.
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
    private static StringBuilder Run(String tagfile, String[] PrimersList, String[] PrimersName, String[] PrimersOri, boolean isprobe, boolean ispattern, boolean iscircle, boolean seqextract, int flanks, int Err3end, int minlen, int maxlen, boolean FRpairs, boolean pcr_predict, boolean alignment, boolean ShowOnlyAmplicons, boolean PCRmatch_alignment, boolean PrimerStatistic) {
        try {
            System.out.println("Running...");
            System.out.println("\nTarget file name: " + tagfile);

            InSilicoPCR s2 = new InSilicoPCR();
            ReadingSequencesFiles rf = new ReadingSequencesFiles(Paths.get(tagfile));
            if (rf.getNseq() == 0) {
                System.out.println("No '>' header line found in " + tagfile
                        + " - the target must be in FASTA format.");
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
            s2.SetSequenceExtractFlanks(flanks);
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

            // The engine's own buffer is handed back rather than copied into a second one: on a
            // genome-sized report that copy doubled the peak heap. Safe only while s2 is local to
            // this method and dies with it - do not hoist the InSilicoPCR instance out.
            StringBuilder sr = s2.getResult();
            sr.append(st);
            return sr;

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
            return null; // load failure -> signal "not processed" to the caller
        }
    }

    // Write text to a file as UTF-8, creating parent directories and truncating any existing file.
    private static void writeFile(String outfile, CharSequence text) throws IOException {
        Path path = Paths.get(outfile);
        Path parent = path.getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        try (BufferedWriter w = Files.newBufferedWriter(path, StandardCharsets.UTF_8)) {
            writeChunked(w, text);
        }
    }

    // Copies the text out in fixed slices. Writer.append(CharSequence) would call toString() on
    // the whole thing, which is the copy this is meant to avoid - and which also caps the report
    // at ~2^31 characters no matter how large the heap is.
    private static void writeChunked(BufferedWriter w, CharSequence text) throws IOException {
        final int step = 1 << 16;
        final int len = text.length();
        char[] buf = new char[Math.min(step, Math.max(len, 1))];
        for (int off = 0; off < len; off += step) {
            int n = Math.min(step, len - off);
            if (text instanceof StringBuilder sb) {
                sb.getChars(off, off + n, buf, 0);
            } else if (text instanceof String str) {
                str.getChars(off, off + n, buf, 0);
            } else {
                for (int i = 0; i < n; i++) {
                    buf[i] = text.charAt(off + i);
                }
            }
            w.write(buf, 0, n);
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
