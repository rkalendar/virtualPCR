
import java.util.ArrayList;
import java.util.List;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

public class virtualpcr {

    public static void main(String[] args) {
        if (args.length > 0) {
            String tagfile = "";
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
            boolean FRpairs = false;              //SetLookR_Fprimes
            boolean pcr_predict = true;           //SetShowPCRProducts
            boolean CalculatePCRproduct = false;  //SetShowPCRproductCalculation
            boolean alignment = true;             //SetShowPrimerAlignment
            boolean ShowOnlyAmplicons = false;    //SetShowOnlyAmplicons
            boolean PCRmatch_alignment = true;    //SetShowPrimerAlignmentPCRproduct
            boolean CTconversion = false;           //CTconversion=false/true

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            try (BufferedReader br = new BufferedReader(new FileReader(infile))) {
                String line;
                while ((line = br.readLine()) != null) {
                    System.out.println(line);
                    line = line.toLowerCase();
                    if (line.contains("targets_path=")) {
                        tagfile = line.substring(13);
                    }
                    if (line.contains("primers_path=")) {
                        primersfile = line.substring(13);
                    }
                    if (line.contains("type=probe")) {
                        isprobe = true;
                    }
                    if (line.contains("molecular=circle")) {
                        iscircle = true;
                    }
                    if (line.contains("linkedsearch=true")) {
                        ispattern = true;
                    }
                    if (line.contains("showprimeralignmentpcrproduct=true")) {
                        PCRmatch_alignment = true;
                    }
                    if (line.contains("sequenceextract=true")) {
                        seqextract = true;
                    }
                    if (line.contains("frpairs=true")) {
                        FRpairs = true;
                    }
                    if (line.contains("showprimeralignment=false")) {
                        alignment = false;
                    }
                    if (line.contains("showpcrproducts=false")) {
                        pcr_predict = false;
                    }
                    if (line.contains("showpcrproductcalculation=true")) {
                        CalculatePCRproduct = true;
                    }
                    if (line.contains("showonlyamplicons=true")) {
                        ShowOnlyAmplicons = true;
                    }
                    if (line.contains("ctconversion=true")) {
                        CTconversion = true;
                    }
                    if (line.contains("number3errors=")) {
                        int h = StrToInt(line.substring(14));
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

            String[] PrimersList = new String[1];
            String[] PrimersName = new String[1];
            String[] PrimersOri = new String[1];
            int n = 1;
            int t = 0;// t=1 Fasta; t=2 TAB

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
                                ln.add(line.trim());
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
                                    if (s.length() > 4) {
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
                    if (folder.isDirectory()) {
                        File[] files = folder.listFiles();
                        int k = -1;
                        String[] filelist = new String[files.length];
                        for (File file : files) {
                            if (file.isFile()) {
                                filelist[++k] = file.getAbsolutePath();
                            }
                        }
                        for (String nfile : filelist) {
                            try {
                                Run(nfile, PrimersList, PrimersName, PrimersOri, isprobe, ispattern, iscircle, seqextract, Err3end, minlen, maxlen, FRpairs, pcr_predict, CalculatePCRproduct, alignment, ShowOnlyAmplicons, PCRmatch_alignment, CTconversion);
                            } catch (Exception e) {
                                System.out.println("Failed to open file: " + nfile);
                            }
                        }

                    } else {
                        Run(tagfile, PrimersList, PrimersName, PrimersOri, isprobe, ispattern, iscircle, seqextract, Err3end, minlen, maxlen, FRpairs, pcr_predict, CalculatePCRproduct, alignment, ShowOnlyAmplicons, PCRmatch_alignment, CTconversion);
                    }
                } else {
                    System.out.println("\nFailed to open file: " + folder);
                }
            } else {
                System.out.println("\nFailed to open primer's file: " + primersfile);
            }
        } else {
            System.out.println("virtualPCR (2024) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/virtualPCR\n");
            System.out.println("Basic usage:");
            System.out.println("java -jar \\virtualPCR\\dist\\virtualPCR.jar \\virtualPCR\\task\\config.file");
        }
    }

    private static void Run(String tagfile, String[] PrimersList, String[] PrimersName, String[] PrimersOri, boolean isprobe, boolean ispattern, boolean iscircle, boolean seqextract, int Err3end, int minlen, int maxlen, boolean FRpairs, boolean pcr_predict, boolean CalculatePCRproduct, boolean alignment, boolean ShowOnlyAmplicons, boolean PCRmatch_alignment, boolean CTconversion) {
        String outfile = tagfile + ".out";
        try {

            StringBuilder sr = new StringBuilder(100000);
            //InSilicoPCR2 s2 = new InSilicoPCR2();
            InSilicoPCR3 s2 = new InSilicoPCR3();

            byte[] binaryArray = Files.readAllBytes(Paths.get(tagfile));
            ReadingSequencesFiles rf = new ReadingSequencesFiles(binaryArray);
            if (rf.getNseq() == 0) {
                System.out.println("There is no sequence(s).");
                System.out.println("File format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
                System.out.println(">seq1\nactacatactacatcactctctctccgcacag\n");
                return;
            }
            System.out.println("\nTarget file name: " + tagfile);
            System.out.println("Target sequence length = " + rf.getLength() + " nt");
            System.out.println("Running...");
            //s2.SetSequences(rf.getSequences(), rf.getNames());
            s2.SetSequences(rf.getSequences(), rf.getNames(), CTconversion);

            s2.SetShowPrimerAlignment(true);
            s2.SetShowPrimerAlignmentPCRproduct(true);
            s2.SetShowOnlyAmplicons(false);
            //s2.SetCTbisulfate(CTconversion);
            s2.SetLookR_Fprimes(FRpairs);
            s2.SetShowPCRProducts(pcr_predict);
            //s2.SetShowPCRproductCalculation(CalculatePCRproduct);
            s2.SetShowPrimerAlignment(alignment);
            s2.SetShowOnlyAmplicons(ShowOnlyAmplicons);
            s2.SetShowPrimerAlignmentPCRproduct(PCRmatch_alignment);
            s2.SetProductMaxLength(maxlen);
            s2.SetProductMinLength(minlen);
            s2.SetPrimers(PrimersList, PrimersName, PrimersOri, isprobe, ispattern, iscircle, seqextract, Err3end);
            long startTime = System.nanoTime();
            s2.Run();
            long duration = (System.nanoTime() - startTime) / 1000000000;
            sr.append(s2.getResult());
            System.out.println(sr);
            System.out.println("Time taken: " + duration + " seconds\n\n");
            try (FileWriter fileWriter = new FileWriter(outfile)) {
                System.out.println("Saving report file: " + outfile);
                fileWriter.write(sr.toString());
                fileWriter.write("Time taken: " + duration + " seconds\n\n");
            }

        } catch (IOException e) {
            System.out.println("Incorrect file name.\n");
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
