import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class PrimerLoader {

    /** Контейнер с загруженными праймерами.
     *  Массивы 1-индексированные: индекс 0 — пустая строка (как в исходном коде). */
    public static class Primers {
        public final String[] list;       // обработанные последовательности (после dna.DNA)
        public final String[] names;      // ID праймеров
        public final String[] originals;  // исходные строки последовательностей
        public final int count;           // количество реально загруженных праймеров

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
     * Загружает праймеры из файла. Формат определяется автоматически:
     *   - FASTA  (есть строка с '>')
     *   - TAB    : ID<TAB>seq1[<TAB>seq2 ...]
     *   - SPACE  : ID seq1 [seq2 ...]
     *
     * В колоночных форматах первая колонка — ID, все остальные — последовательности.
     * Последовательности длиной <= 11 после dna.DNA(...) пропускаются.
     * @param primersfile
     */
    public static Primers loadPrimers(String primersfile) {

        // --- 1. Прочитать файл ---
        List<String> lines;
        try {
            lines = Files.readAllLines(Paths.get(primersfile));
        } catch (IOException e) {
            System.out.println("\nFailed to open primer's file: " + primersfile);
            return Primers.empty();
        }

        // --- 2. Определить формат: 0=space, 1=FASTA, 2=TAB ---
        int format = 0;
        for (String line : lines) {
            if (line.trim().startsWith(">"))  { format = 1; break; } // a stray '>' inside a column must not flip the whole file to FASTA
            if (line.contains("\t")) { format = 2; break; }
        }

        // --- 3. Распарсить ---
        List<String> names     = new ArrayList<>();
        List<String> originals = new ArrayList<>();
        List<String> processed = new ArrayList<>();

        String currentName = null;
        StringBuilder currentSeq = new StringBuilder();

        for (String line : lines) {
            if (line.isEmpty()) continue;

            if (format == 1) {                  // FASTA
                String trimmed = line.trim();
                if (trimmed.startsWith(">")) {
                    flushFasta(currentName, currentSeq, names, originals, processed);
                    currentName = trimmed.substring(1).trim();
                    currentSeq.setLength(0);
                } else {
                    currentSeq.append(trimmed);
                }
            } else {                            // колоночный формат
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

        // --- 4. Сформировать 1-индексированные массивы ---
        int n = processed.size();
        if (n == 0) return Primers.empty();

        String[] PrimersList = new String[n + 1];
        String[] PrimersName = new String[n + 1];
        String[] PrimersOri  = new String[n + 1];
        PrimersList[0] = "";
        PrimersName[0] = "";
        PrimersOri[0]  = "";
        for (int i = 0; i < n; i++) {
            PrimersList[i + 1] = processed.get(i);
            PrimersName[i + 1] = names.get(i);
            PrimersOri[i + 1]  = originals.get(i);
        }
        return new Primers(PrimersList, PrimersName, PrimersOri);
    }

    // ---- helpers ----

    private static void flushFasta(String name, StringBuilder seq,
                                   List<String> names, List<String> originals,
                                   List<String> processed) {
        if (name == null || seq.length() == 0) return;
        addPrimer(name, seq.toString(), names, originals, processed);
    }

    private static void addPrimer(String id, String rawSeq,
                                  List<String> names, List<String> originals,
                                  List<String> processed) {
        String s = dna.DNA(rawSeq.toLowerCase());   // прямой статический вызов
        if (s.length() > 11) {
            names.add(id);
            originals.add(rawSeq);
            processed.add(s);
        }
    }
}