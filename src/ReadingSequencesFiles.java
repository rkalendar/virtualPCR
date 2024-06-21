
import java.util.Arrays;

public final class ReadingSequencesFiles {

    public ReadingSequencesFiles(byte[] s) {
        source = s;
        ReadingSequences();
    }

    public ReadingSequencesFiles(String s) {
        source = s.getBytes();
        ReadingSequences();
    }

    public String[] getSequences() {
        if (ns == 0) {
            return null;
        }
        return sequence;
    }

    public String[] getNames() {
        if (ns == 0) {
            return null;
        }
        return name_seq;
    }

    public int getNseq() {
        return ns;
    }

    private void ReadingSequences() {
        if (source == null) {
            return;
        }
        dnay = new double[128];

        int l = source.length;
        int s = 0; // total length
        ns = 0;    // amount fasta sequences

        for (int i = 0; i < l; i++) {
            if (source[i] < 9) {
                ns = 0;
                break;
            }
            if (tables.cdn[source[i]] > 0) {
                s++;
            }
            if (source[i] == 62) {
                ns++;
            }
        }
        if (s == 0 || ns == 0) {
            return;
        }
        name_seq = new String[ns];
        sequence = new String[ns];
        int n = -1;
        int t = 0;

        for (int i = 0; i < l; i++) {
            if (source[i] == 62) {
                if (t > 0) {

                    byte[] d = Arrays.copyOfRange(source, t, i - 1); //public static short[] copyOfRange(short[] original, int from, int to)
                    int x = 0;
                    for (int j = 0; j < d.length; j++) {
                        if (tables.cdn[d[j]] > 0) {
                            d[x] = tables.cdn[d[j]];
                            x++;
                        }
                    }
                    sequence[n] = new String(Arrays.copyOfRange(d, 0, x));
                    lSeqs = lSeqs + x;
                }
                n++;
                for (int j = i + 1; j < l; j++) {
                    if (source[j] == 10 || source[j] == 13) {
                        name_seq[n] = new String(source, i + 1, j - i).trim();
                        i = j;
                        t = j + 1;
                        break;
                    }
                }
            }
        }
        byte[] d = Arrays.copyOfRange(source, t, l);
        int x = 0;
        for (int j = 0; j < d.length; j++) {
            if (tables.cdn[d[j]] > 0) {
                d[x] = tables.cdn[d[j]];
                x++;
            }
        }
        lSeqs = lSeqs + x;
        sequence[n] = new String(Arrays.copyOfRange(d, 0, x));
    }

    public void SetFolder(String p) {
        if (!p.isEmpty()) {
            fld = p + "__";
        }
    }

    public double getA() {
        return (dnay[97] + (dnay[109] + dnay[114] + dnay[119]) / 2 + (dnay[118] + dnay[104] + dnay[100]) / 3 + dnay[110] / 4);
    }

    public double getT() {
        return (dnay[116] + dnay[117] + (dnay[121] + dnay[107] + dnay[119]) / 2 + (dnay[98] + dnay[104] + dnay[100]) / 3 + dnay[110] / 4);
    }

    public double getC() {
        return (dnay[99] + (dnay[109] + dnay[115] + dnay[121]) / 2 + (dnay[98] + dnay[118] + dnay[104]) / 3 + dnay[110] / 4);
    }

    public double getG() {
        return (dnay[103] + dnay[105] + ((dnay[115] + dnay[114] + dnay[107]) / 2) + ((dnay[98] + dnay[118] + dnay[100]) / 3) + dnay[110] / 4);
    }

    public double getN() {
        return dnay[110];
    }

    public double getR() {
        return (dnay[103] + dnay[105] + dnay[97] + dnay[114] + ((dnay[115] + dnay[107] + dnay[109] + dnay[119] + dnay[110]) / 2) + ((dnay[98] + dnay[118] + dnay[118] + dnay[104] + dnay[100] + dnay[100]) / 3));
    }

    public double getY() {
        return (dnay[99] + dnay[116] + dnay[117] + dnay[121] + (dnay[109] + dnay[115] + dnay[110] + dnay[107] + dnay[119]) / 2 + (dnay[98] + dnay[118] + dnay[104] + dnay[98] + dnay[104] + dnay[100]) / 3);
    }

    public int getLength() {
        return (int) lSeqs;
    }

    public double getCG() {
        return lSeqs < 1 ? 0 : (100 * (dnay[103] + dnay[105] + dnay[99] + dnay[115] + ((dnay[114] + dnay[107] + dnay[109] + dnay[121] + dnay[110]) / 2) + ((dnay[98] + dnay[118] + dnay[100] + dnay[98] + dnay[118] + dnay[104]) / 3))) / lSeqs;
    }
    private String fld = "";
    private double[] dnay;
    private String[] name_seq;
    private String[] sequence;
    private byte[] source = null;
    private double lSeqs = 0;
    private int ns = 0;
}
