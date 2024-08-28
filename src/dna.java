
public final class dna {

    final public static int maxdna = 1025; // =2^12 or 4^6=4096 (+1); 3^12=531,441; 4^12=16,777,216
    final public static int mindna = 1025; // 2^10 - 4^5 

    public static String AntisenseDNA(String source) {
        byte[] cdna = new byte[128];
        cdna[65] = 84;  //A
        cdna[66] = 86;  //B
        cdna[67] = 71;  //C
        cdna[68] = 72;  //D
        cdna[71] = 67;  //G
        cdna[72] = 68;  //H
        cdna[73] = 71;  //I
        cdna[75] = 77;  //K
        cdna[77] = 75;  //M
        cdna[78] = 78;  //N
        cdna[82] = 89;  //R
        cdna[83] = 83;  //S
        cdna[84] = 65;  //T
        cdna[85] = 65;  //U
        cdna[86] = 66;  //V
        cdna[87] = 87;  //W
        cdna[89] = 82;  //Y

        cdna[97] = 116;//  t <- a
        cdna[98] = 118;//  v <- b
        cdna[99] = 103;//  g <- c
        cdna[100] = 104;// h <- d
        cdna[103] = 99; // c <- g
        cdna[104] = 100;// d <- h
        cdna[105] = 99; // i <- g
        cdna[107] = 109;// m <- k
        cdna[109] = 107;// k <- m
        cdna[110] = 110;// n <- n
        cdna[114] = 121;// y <- r
        cdna[115] = 115;// s <- s
        cdna[116] = 97; // a <- t
        cdna[117] = 97; // a <- u
        cdna[118] = 98; // b <- v
        cdna[119] = 119;// w <- w
        cdna[121] = 114;// r <- y

        byte[] b = source.getBytes();
        for (int i = 0; i < source.length(); i++) {
            if (cdna[b[i]] > 0) {
                if (b[i] < 97) {
                    b[i] = (byte) (b[i] + 32);
                }
                b[i] = cdna[b[i]];
            }
        }
        return new String(b);
    }

    public static String ReverseSeq(String source) {
        return new StringBuilder(source).reverse().toString();
    }

    public static String ReverseSeq(char[] source) {
        StringBuilder s = new StringBuilder();
        return s.append(source).reverse().toString();
    }

    public static String ComplementDNA(String source) {
        byte[] cdna = new byte[128];
        cdna[65] = 84;  //A
        cdna[66] = 86;  //B
        cdna[67] = 71;  //C
        cdna[68] = 72;  //D
        cdna[71] = 67;  //G
        cdna[72] = 68;  //H
        cdna[73] = 71;  //I
        cdna[75] = 77;  //K
        cdna[77] = 75;  //M
        cdna[78] = 78;  //N
        cdna[82] = 89;  //R
        cdna[83] = 83;  //S
        cdna[84] = 65;  //T
        cdna[85] = 65;  //U
        cdna[86] = 66;  //V
        cdna[87] = 87;  //W
        cdna[89] = 82;  //Y

        cdna[97] = 116;//  t <- a
        cdna[98] = 118;//  v <- b
        cdna[99] = 103;//  g <- c
        cdna[100] = 104;// h <- d
        cdna[103] = 99; // c <- g
        cdna[104] = 100;// d <- h
        cdna[105] = 99; // i <- g
        cdna[107] = 109;// m <- k
        cdna[109] = 107;// k <- m
        cdna[110] = 110;// n <- n
        cdna[114] = 121;// y <- r
        cdna[115] = 115;// s <- s
        cdna[116] = 97; // a <- t
        cdna[117] = 97; // a <- u
        cdna[118] = 98; // b <- v
        cdna[119] = 119;// w <- w
        cdna[121] = 114;// r <- y

        byte[] b = source.getBytes();
        int l = source.length();
        int n = l / 2;
        for (int i = 0; i < n; i++) {
            if (cdna[b[i]] > 0) {
                byte t = cdna[b[l - i - 1]];
                b[l - i - 1] = cdna[b[i]];
                b[i] = t;
            }
        }
        if ((l % 2) == 1) {
            if (cdna[b[n]] > 0) {
                b[n] = cdna[b[n]];
            }
        }
        return new String(b);
    }

    public static int DNAdegenerateN(String source) {
        if (source == null || source.isEmpty()) {
            return 0;
        }
        int l = source.length();
        if (l < 1) {
            return 0;
        }
        int n = 1;
        byte[] f = source.getBytes();
        for (int i = 0; i < l; i++) {
            n = n * tables.t2[f[i]][0];
            if (n > maxdna) {
                break;
            }
        }
        return n;
    }

    public static String DNActBisulfite(String source) {
        if (source == null || source.isEmpty()) {
            return source;
        }
        String s = source.replaceAll("cg", "11");
        s = s.replaceAll("c", "t");
        s = s.replaceAll("11", "cg");
        return s;

        /*
seq = Replace(s, "cg", "11")
seq = Replace(seq, "c", "t")
DNA_bisulfite = Replace(seq, "11", "cg")       

           s = "3'-" + dna.ReverseSeq(dna.DNActBisulfite(dna.ComplementDNA(Seq[i])));
       
         */
    }

    public static String DNActBisulfite2(String source) {
        if (source == null || source.isEmpty()) {
            return source;
        }
        return DNActBisulfite(ComplementDNA(source));
    }

    public static byte[][] ConvertDNAb(String source) {
        if (source == null || source.isEmpty()) {
            return null;
        }
        int l = source.length();
        if (l < 1) {
            return null;
        }

        int n = 1;
        int i = 0;
        int j = 0;

        byte[] f = source.getBytes();
        byte[][] r = new byte[1][];

        for (i = 0; i < l; i++) {
            n = n * tables.t2[f[i]][0];
            if (n > maxdna) {
                break;
            }
        }
        if (n == 1 || n > maxdna) {
            r[0] = f;
            for (i = 0; i < l; i++) {
                r[0][i] = tables.t2[f[i]][1];
            }
            return (r);
        }

        r = new byte[n][];
        r[0] = f.clone();

        for (i = 0; i < l; i++) {
            r[0][i] = tables.t2[f[i]][1];
        }

        int k = 1;
        for (i = 0; i < l; i++) {
            int h = tables.t2[f[i]][0];
            if (h > 1) {
                for (j = 0; j < h; j++) {
                    int x = 0;
                    byte z = tables.t2[f[i]][j + 1];
                    for (int p = j * k; p < j * k + k; p++) {
                        r[p] = r[x++].clone();
                        r[p][i] = z;
                    }
                }
                k = k * h;
            }
        }
        return r;
    }

    public static String[] ConvertDNA(String source) {
        if (source == null || source.isEmpty()) {
            return null;
        }
        int l = source.length();
        if (l < 1) {
            return null;
        }

        int n = 1;
        int i = 0;
        int j = 0;

        String[] s = new String[]{source};

        byte[] f = source.getBytes();
        byte[][] r = new byte[1][];

        for (i = 0; i < l; i++) {
            n = n * tables.t3[f[i]][0];
            if (n > maxdna) {
                break;
            }
        }
        if (n == 1 || n > maxdna) {
            r[0] = f;
            for (i = 0; i < l; i++) {
                r[0][i] = tables.t3[f[i]][1];
            }
            return (s);
        }

        r = new byte[n][];
        r[0] = f.clone();

        int k = 1;
        for (i = 0; i < l; i++) {
            int h = tables.t3[f[i]][0];
            if (h > 1) {
                for (j = 0; j < h; j++) {
                    int x = 0;
                    byte z = tables.t3[f[i]][j + 1];
                    for (int p = j * k; p < j * k + k; p++) {
                        r[p] = r[x++].clone();
                        r[p][i] = z;
                    }
                }
                k = k * h;
            }
        }

        s = new String[n];
        for (int h = 0; h < n; h++) {
            s[h] = new String(r[h]);
        }
        return s;
    }

    public static String DNA(byte[] b) {
        int l = b.length;
        byte[] cdn = new byte[128];
        cdn[65] = 2;  // A
        cdn[66] = 2;  // B
        cdn[67] = 2;  // C
        cdn[68] = 2;  // D
        cdn[71] = 2;  // G
        cdn[72] = 2;  // H
        cdn[73] = 2;  // I
        cdn[75] = 2;  // K
        cdn[77] = 2;  // M
        cdn[78] = 2;  // N
        cdn[82] = 2;  // R
        cdn[83] = 2;  // S
        cdn[84] = 2;  // T
        cdn[85] = 2;  // U
        cdn[86] = 2;  // V
        cdn[87] = 2;  // W
        cdn[89] = 2;  // Y        
        cdn[97] = 1;   // a
        cdn[98] = 1;   // b
        cdn[99] = 1;   // c
        cdn[100] = 1;  // d
        cdn[103] = 1;  // g
        cdn[104] = 1;  // h
        cdn[105] = 1;  // i
        cdn[107] = 1;  // k
        cdn[109] = 1;  // m
        cdn[110] = 1;  // n
        cdn[114] = 1;  // r
        cdn[115] = 1;  // s
        cdn[116] = 1;  // t
        cdn[117] = 1;  // u
        cdn[118] = 1;  // v
        cdn[119] = 1;  // w
        cdn[121] = 1;  // y  

        int n = -1;
        for (int i = 0; i < l; i++) {
            if (b[i] < 122 && b[i] > 64) {
                int t = cdn[b[i]];
                if (t > 0) {
                    b[++n] = (t == 2) ? (byte) (b[i] + 32) : b[i];
                }
            }
        }
        return (n > -1) ? (new String(b, 0, n + 1)) : "";
    }

    public static String DNA(String source) {
        if (source == null || source.isEmpty()) {
            return "";
        }
        byte[] cdn = new byte[128];
        cdn[65] = 2;  // A
        cdn[66] = 2;  // B
        cdn[67] = 2;  // C
        cdn[68] = 2;  // D
        cdn[71] = 2;  // G
        cdn[72] = 2;  // H
        cdn[73] = 2;  // I
        cdn[75] = 2;  // K
        cdn[77] = 2;  // M
        cdn[78] = 2;  // N
        cdn[82] = 2;  // R
        cdn[83] = 2;  // S
        cdn[84] = 2;  // T
        cdn[85] = 2;  // U
        cdn[86] = 2;  // V
        cdn[87] = 2;  // W
        cdn[89] = 2;  // Y        
        cdn[97] = 1;   // a
        cdn[98] = 1;   // b
        cdn[99] = 1;   // c
        cdn[100] = 1;  // d
        cdn[103] = 1;  // g
        cdn[104] = 1;  // h
        cdn[105] = 1;  // i
        cdn[107] = 1;  // k
        cdn[109] = 1;  // m
        cdn[110] = 1;  // n
        cdn[114] = 1;  // r
        cdn[115] = 1;  // s
        cdn[116] = 1;  // t
        cdn[117] = 1;  // u
        cdn[118] = 1;  // v
        cdn[119] = 1;  // w
        cdn[121] = 1;  // y  

        int l = source.length();
        byte[] b = source.getBytes();
        int n = -1;
        for (int i = 0; i < l; i++) {
            if (b[i] < 122 && b[i] > 64) {
                int t = cdn[b[i]];
                if (t > 0) {
                    b[++n] = (t == 2) ? (byte) (b[i] + 32) : b[i];
                }
            }
        }
        return (n > -1) ? (new String(b, 0, n + 1)) : "";
    }
}
