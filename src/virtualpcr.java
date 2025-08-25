import java.util.ArrayList;
import java.util.List;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.InvalidPathException;
import java.nio.file.Path;
import java.nio.file.Paths;
import static java.nio.file.StandardOpenOption.APPEND;
import static java.nio.file.StandardOpenOption.CREATE;

public class virtualpcr {

    public static void main(String[] args) {
        if (args.length > 0) {
            String tagfile = "";
            String outfile = "";
            String primersfile = "";
            String infile = args[0];
            for (String arg : args) {
                if (args.length > 0) {
                    infile = arg;
                    break;
                }
            }

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
                        if (value.length() > 1) {
                            outfile = value;
                        }
                    }

                    if (line.contains("targets_path=")) {
                        tagfile = cline.substring(13).trim();
                        System.out.println(tagfile);
                    }
                    if (line.contains("primers_path=")) {
                        primersfile = cline.substring(13).trim();
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
                        int h = StrToInt(line.substring(14));
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
                        int h = StrToInt(line.substring(6));
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
                        int h = StrToInt(line.substring(6));
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
            }
            System.out.println("\n");

            String[] PrimersList = new String[1];
            String[] PrimersName = new String[1];
            String[] PrimersOri = new String[1];
            int n = 1;
            int t = 0;//t=0 Space; t=1 Fasta; t=2 TAB 

            try {
                List<String> lines = Files.readAllLines(Paths.get(primersfile));
                for (String line : lines) {
                    if (line.contains(">")) {
                        t = 1;
                        break;
                    }
                    if (line.contains("\t")) {
                        t = 2;
                        break;
                    }
                }
            } catch (IOException e) {
                System.out.println("\nFailed to open primer's file: " + primersfile);
                return;
            }

            try {
                List<String> ln = new ArrayList<>();
                List<String> lo = new ArrayList<>();
                List<String> ls = new ArrayList<>();
                List<String> lines = Files.readAllLines(Paths.get(primersfile));
                for (String line : lines) {
                    if (!line.isEmpty()) {
                        System.out.println(line.trim());
                        if (t == 1) {
                            if (line.contains(">")) {
                                n++;
                                ln.add(line.substring(1).trim());
                            } else {
                                lo.add(line.trim());
                                ls.add(dna.DNA(line.toLowerCase()));
                            }
                        }
                        if (t == 2) {
                            String[] a = line.split("[ \t]+");
                            if (a.length > 1) {
                                for (int i = 1; i < a.length; i++) {
                                    String s = dna.DNA(a[i].toLowerCase());
                                    if (s.length() > 11) {
                                        ln.add(a[0]);
                                        lo.add(a[i]);
                                        ls.add(s);
                                        n++;
                                    }
                                }
                            }
                        }
                        if (t == 0) {
                            String[] a = line.split("[ ]+");
                            if (a.length > 1) {
                                for (int i = 1; i < a.length; i++) {
                                    String s = dna.DNA(a[i].toLowerCase());
                                    if (s.length() > 11) {
                                        ln.add(a[0]);
                                        lo.add(a[i]);
                                        ls.add(s);
                                        n++;
                                    }
                                }
                            }
                        }
                    }
                }

                if (n > 1) {
                    PrimersList = new String[n];
                    PrimersName = new String[n];
                    PrimersOri = new String[n];
                    PrimersList[0] = "";
                    PrimersName[0] = "";
                    PrimersOri[0] = "";
                    for (int i = 0; i < n - 1; i++) {
                        PrimersList[i + 1] = ls.get(i);
                        PrimersName[i + 1] = ln.get(i);
                        PrimersOri[i + 1] = lo.get(i);
                    }
                }
            } catch (IOException e) {
                return;
            }

            if (n > 1) {
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
                        int k = -1;
                        String[] filelist = new String[files.length];
                        for (File file : files) {
                            if (file.isFile()) {
                                filelist[++k] = file.getAbsolutePath();
                            }
                        }

                        try (FileWriter fileWriter1 = new FileWriter(outfile)) {
                            for (String nfile : filelist) {
                                try {
                                    StringBuilder sr = Run(nfile, outfile, PrimersList, PrimersName, PrimersOri, isprobe, ispattern, iscircle, seqextract, Err3end, minlen, maxlen, FRpairs, pcr_predict, alignment, ShowOnlyAmplicons, PCRmatch_alignment, CTconversion, PrimerStatistic);
                                    fileWriter1.write(sr.toString());
                                    fileWriter1.write("\n\n");
                                } catch (IOException e) {
                                    System.out.println("Failed to open file: " + nfile);
                                }
                            }
                        } catch (IOException e) {
                            System.out.println("Incorrect file name.\n");
                        }

                    } else {
                        StringBuilder sr = Run(tagfile, outfile, PrimersList, PrimersName, PrimersOri, isprobe, ispattern, iscircle, seqextract, Err3end, minlen, maxlen, FRpairs, pcr_predict, alignment, ShowOnlyAmplicons, PCRmatch_alignment, CTconversion, PrimerStatistic);
                    }
                } else {
                    System.out.println("\nFailed to open file: " + folder);
                }
            } else {
                System.out.println("\nFailed to open primer's file: " + primersfile);
            }
        } else {
            System.out.println("virtualPCR (2024-2025) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/virtualPCR\n");
            System.out.println("Basic usage:");
            System.out.println("java -jar \\virtualPCR\\dist\\virtualPCR.jar \\virtualPCR\\task\\config.file");
        }
    }

    private static StringBuilder Run(String tagfile, String outfile, String[] PrimersList, String[] PrimersName, String[] PrimersOri, boolean isprobe, boolean ispattern, boolean iscircle, boolean seqextract, int Err3end, int minlen, int maxlen, boolean FRpairs, boolean pcr_predict, boolean alignment, boolean ShowOnlyAmplicons, boolean PCRmatch_alignment, boolean CTconversion, boolean PrimerStatistic) {
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
                return sr;
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
            appendLine(outfile, sr.toString());
            System.out.println(sr);
            System.out.println("Saved report: " + outfile);

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
        }
        return sr;
    }

    public static void appendLine(String outfile, CharSequence text) throws IOException {
        Path path = Paths.get(outfile);
        Path parent = path.getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        try (BufferedWriter w = Files.newBufferedWriter(path, StandardCharsets.UTF_8, CREATE, APPEND)) {
            w.append(text);
            w.append(System.lineSeparator());
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
        return (Integer.parseInt(r.toString()));
    }
}
