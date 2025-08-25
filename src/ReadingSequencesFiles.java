import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

public final class ReadingSequencesFiles {

    private String[] name_seq;
    private String[] sequence;
    private int lSeqs = 0;
    private int ns = 0;
    private byte[] dnl;

    public ReadingSequencesFiles(Path fastaPath) throws IOException {
        LoadTable();
        readFastaStream(Files.newBufferedReader(fastaPath, StandardCharsets.US_ASCII), /*normalize*/ true);
    }

    public static ReadingSequencesFiles readMasking(Path fastaPath) throws IOException {
        ReadingSequencesFiles r = new ReadingSequencesFiles();
        r.readFastaStream(Files.newBufferedReader(fastaPath, StandardCharsets.US_ASCII), /*normalize*/ false);
        return r;
    }

    public ReadingSequencesFiles(byte[] s) {
        LoadTable();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(s), StandardCharsets.US_ASCII))) { //
            readFastaStream(br, true);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public ReadingSequencesFiles(String s) {
        LoadTable();
        try (BufferedReader br = new BufferedReader(new StringReader(s))) {
            readFastaStream(br, true);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public ReadingSequencesFiles() {
        LoadTable();
    }

    public String[] getSequences() {
        return ns == 0 ? null : sequence;
    }

    public String[] getNames() {
        return ns == 0 ? null : name_seq;
    }

    public int getNseq() {
        return ns;
    }

    public int getLength() {
        return lSeqs;
    }

    private void LoadTable() {
        dnl = new byte[128]; // ASCII
        dnl[65] = 97;  // A->a
        dnl[66] = 98;   // B->b
        dnl[67] = 99;   // C->c
        dnl[68] = 100;  // D->d
        dnl[71] = 103;  // G->g
        dnl[72] = 104;  // H->h
        dnl[73] = 99;   // I->c  
        dnl[75] = 107;  // K->k
        dnl[77] = 109;  // M->m
        dnl[78] = 110;  // N->n
        dnl[82] = 114;  // R->r
        dnl[83] = 115;  // S->s
        dnl[84] = 116;  // T->t
        dnl[85] = 116;  // U->t
        dnl[86] = 118;  // V->v
        dnl[87] = 119;  // W->w
        dnl[89] = 121;  // Y->y

        dnl[97] = 97;   // a
        dnl[98] = 98;   // b
        dnl[99] = 99;   // c
        dnl[100] = 100;  // d
        dnl[103] = 103;  // g
        dnl[104] = 104;  // h
        dnl[105] = 99;   // i -> c
        dnl[107] = 107;  // k
        dnl[109] = 109;  // m
        dnl[110] = 110;  // n
        dnl[114] = 114;  // r
        dnl[115] = 115;  // s
        dnl[116] = 116;  // t
        dnl[117] = 116;  // u -> t
        dnl[118] = 118;  // v
        dnl[119] = 119;  // w
        dnl[121] = 121;  // y
    }

    private void readFastaStream(BufferedReader br, boolean normalize) throws IOException {
        List<String> names = new ArrayList<>();
        List<String> seqs = new ArrayList<>();

        String line;
        String currentName = null;
        StringBuilder currentSeq = new StringBuilder(1 << 20); 

        while ((line = br.readLine()) != null) {
            if (line.isEmpty()) {
                continue;
            }

            if (line.charAt(0) == '>') {
                if (currentName != null) {
                    seqs.add(currentSeq.toString());
                    currentSeq.setLength(0);
                }
                currentName = line.substring(1).trim();
                names.add(currentName);
                continue;
            }

            final int len = line.length();
            for (int i = 0; i < len; i++) {
                char ch = line.charAt(i);
                if (ch < 128 && dnl[ch] > 0) {
                    if (normalize) {
                        currentSeq.append((char) dnl[ch]);
                    } else {
                        currentSeq.append(ch);
                    }
                    lSeqs++;
                }

            }
        }

        if (currentName != null) {
            seqs.add(currentSeq.toString());
        }

        ns = seqs.size();
        if (ns == 0) {
            name_seq = null;
            sequence = null;
            lSeqs = 0;
            return;
        }

        name_seq = names.toArray(String[]::new);
        sequence = seqs.toArray(String[]::new);
    }

    public void ReadingMaskSequencesFiles(byte[] source) {
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(source), StandardCharsets.US_ASCII))) {
            this.name_seq = null;
            this.sequence = null;
            this.ns = 0;
            this.lSeqs = 0;
            readFastaStream(br, false);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

}
