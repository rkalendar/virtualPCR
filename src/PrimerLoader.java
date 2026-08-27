import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class PrimerLoader {

    /**
     * Container containing loaded primers. Arrays are 1-indexed: index 0
     * — an empty string (as in the source code).
     */
    public static class Primers {

        public final String[] list;       // processed sequences (after dna.DNA)
        public final String[] names;      // ID primers
        public final String[] originals;  // initial sequence strings
        public final int count;           // the number of primers actually loaded

        public Primers(String[] list, String[] names, String[] originals) {
            this.list = list;
            this.names = names;
            this.originals = originals;
            this.count = list.length - 1;
        }

        public static Primers empty() {
            return new Primers(new String[]{""}, new String[]{""}, new String[]{""});
        }
    }

    /**
     * Loads primers from a file. The format is detected automatically: - FASTA
     * (contains a line with '>') - TAB: ID<TAB>seq1[<TAB>seq2 ...] - SPACE: ID seq1
     * [seq2 ...]
     *
     * In columnar formats, the first column is the ID, and all others are the sequences. Sequences of length <= 11 after dna.DNA(...) are skipped. @param primersfile
     */
    // Reports nothing on failure: the caller knows whether the path was blank, the file missing,
    // or the contents unusable, and prints one accurate message for the case that actually applies.
    public static Primers loadPrimers(String primersfile) {

        // --- 1. Прочитать файл ---
        List<String> lines;
        try {
            lines = Files.readAllLines(Paths.get(primersfile));
        } catch (java.nio.charset.MalformedInputException me) {
            // 8-bit primer list (Windows-1251/Latin-1): ISO-8859-1 maps all 256 byte values, so a
            // hand-edited file with non-ASCII identifiers still loads. Tried second, so a genuine
            // UTF-8 file keeps decoding correctly instead of turning into mojibake.
            // MalformedInputException is an IOException, so this catch must come first.
            try {
                lines = Files.readAllLines(Paths.get(primersfile), java.nio.charset.StandardCharsets.ISO_8859_1);
            } catch (IOException e) {
                return Primers.empty();
            }
        } catch (IOException e) {
            return Primers.empty();
        }

        // --- 2. Specify the format: 0=space, 1=FASTA, 2=TAB ---
        int format = 0;
        for (String line : lines) {
            if (line.trim().startsWith(">")) {
                format = 1;
                break;
            } // a stray '>' inside a column must not flip the whole file to FASTA
            if (line.contains("\t")) {
                format = 2;
                break;
            }
        }

        // --- 3. Parse ---
        List<String> names = new ArrayList<>();
        List<String> originals = new ArrayList<>();
        List<String> processed = new ArrayList<>();

        String currentName = null;
        StringBuilder currentSeq = new StringBuilder();

        for (String line : lines) {
            if (line.isEmpty()) {
                continue;
            }

            if (format == 1) {                  // FASTA
                String trimmed = line.trim();
                if (trimmed.startsWith(">")) {
                    flushFasta(currentName, currentSeq, names, originals, processed);
                    currentName = trimmed.substring(1).trim();
                    currentSeq.setLength(0);
                } else {
                    currentSeq.append(trimmed);
                }
            } else {                            // columnar format
                String[] a = (format == 2) ? line.trim().split("[ \t]+") : line.trim().split("[ ]+"); // trim first so leading whitespace doesn't make an empty ID
                if (a.length > 1) {
                    String id = a[0];
                    for (int i = 1; i < a.length; i++) {
                        addPrimer(id, a[i], names, originals, processed);
                    }
                }
            }
        }
        // последний FASTA-блок
        if (format == 1) {
            flushFasta(currentName, currentSeq, names, originals, processed);
        }

        // --- 4. Create 1-indexed arrays ---
        int n = processed.size();
        if (n == 0) {
            return Primers.empty();
        }

        String[] PrimersList = new String[n + 1];
        String[] PrimersName = new String[n + 1];
        String[] PrimersOri = new String[n + 1];
        PrimersList[0] = "";
        PrimersName[0] = "";
        PrimersOri[0] = "";
        for (int i = 0; i < n; i++) {
            PrimersList[i + 1] = processed.get(i);
            PrimersName[i + 1] = names.get(i);
            PrimersOri[i + 1] = originals.get(i);
        }
        return new Primers(PrimersList, PrimersName, PrimersOri);
    }

    // ---- helpers ----
    private static void flushFasta(String name, StringBuilder seq,
            List<String> names, List<String> originals,
            List<String> processed) {
        if (name == null || seq.length() == 0) {
            return;
        }
        addPrimer(name, seq.toString(), names, originals, processed);
    }

    private static void addPrimer(String id, String rawSeq,
            List<String> names, List<String> originals,
            List<String> processed) {
        String s = dna.DNA(rawSeq.toLowerCase());   // direct static call
        if (s.length() > 5) {
            if (s.length() < 11) {
                s = tools.Strings(12 - s.length(), 'n') + s;
            }
            names.add(id);
            originals.add(rawSeq);
            processed.add(s);
        }
    }
}
