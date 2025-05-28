
public final class dna {

    final public static int maxdna = 1025; // =2^12 or 4^6=4096 (+1); 3^12=531,441; 4^12=16,777,216
    final public static int mindna = 1025; // 2^10 - 4^5 

    public static String AntisenseDNA(String source) {
        byte[] cdna = new byte[128];
        cdna[65] = 116;  //A->t
        cdna[66] = 118;  //B->v
        cdna[67] = 103;  //C->g
        cdna[68] = 104;  //D->h
        cdna[71] = 99;   //G->c
        cdna[72] = 100;  //H->d
        cdna[73] = 99;   //I->c
        cdna[75] = 109;  //K->m
        cdna[77] = 107;  //M->k
        cdna[78] = 110;  //N->n
        cdna[82] = 121;  //R->y
        cdna[83] = 115;  //S->s
        cdna[84] = 97;   //T->a
        cdna[85] = 97;   //U->a
        cdna[86] = 98;   //V->b
        cdna[87] = 119;  //W->w
        cdna[89] = 114;  //Y->r

        cdna[97] = 116;//  t <- a
        cdna[98] = 118;//  v <- b
        cdna[99] = 103;//  g <- c
        cdna[100] = 104;// h <- d
        cdna[103] = 99; // c <- g
        cdna[104] = 100;// d <- h
        cdna[105] = 99; // i -> c
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
        cdna[65] = 116;  //A->t
        cdna[66] = 118;  //B->v
        cdna[67] = 103;  //C->g
        cdna[68] = 104;  //D->h
        cdna[71] = 99;   //G->c
        cdna[72] = 100;  //H->d
        cdna[73] = 99;   //I->c
        cdna[75] = 109;  //K->m
        cdna[77] = 107;  //M->k
        cdna[78] = 110;  //N->n
        cdna[82] = 121;  //R->y
        cdna[83] = 115;  //S->s
        cdna[84] = 97;   //T->a
        cdna[85] = 97;   //U->a
        cdna[86] = 98;   //V->b
        cdna[87] = 119;  //W->w
        cdna[89] = 114;  //Y->r

        cdna[97] = 116;//  t <- a
        cdna[98] = 118;//  v <- b
        cdna[99] = 103;//  g <- c
        cdna[100] = 104;// h <- d
        cdna[103] = 99; // c <- g
        cdna[104] = 100;// d <- h
        cdna[105] = 99; // c <- i
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
        cdn[65] = 97;   // A
        cdn[66] = 98;   // B
        cdn[67] = 99;   // C
        cdn[68] = 100;  // D
        cdn[69] = 97;   // dA=E
        cdn[70] = 99;   // dC=F
        cdn[71] = 103;  // G
        cdn[72] = 104;  // H
        cdn[73] = 103;  // I 
        cdn[74] = 103;  // dG=J
        cdn[75] = 107;  // K
        cdn[76] = 116;  // dT=L   
        cdn[77] = 109;  // M
        cdn[78] = 110;  // N
        cdn[82] = 114;  // R
        cdn[83] = 115;  // S
        cdn[84] = 116;  // T
        cdn[85] = 116;  // U
        cdn[86] = 118;  // V
        cdn[87] = 119;  // W
        cdn[89] = 121;  // Y 

        cdn[97] = 97;    // a
        cdn[98] = 98;    // b
        cdn[99] = 99;    // c
        cdn[100] = 100;  // d         
        cdn[101] = 97;   // dA=E
        cdn[102] = 99;   // dC=F
        cdn[103] = 103;  // g
        cdn[104] = 104;  // h  
        cdn[105] = 103;  // i
        cdn[106] = 103;  // dG=J
        cdn[107] = 107;  // k
        cdn[108] = 116;  // dT=L             
        cdn[109] = 109;  // m
        cdn[110] = 110;  // n
        cdn[114] = 114;  // r
        cdn[115] = 115;  // s
        cdn[116] = 116;  // t 
        cdn[117] = 116;  // u
        cdn[118] = 118;  // v
        cdn[119] = 119;  // w
        cdn[121] = 121;  // y  

        int n = -1;
        for (int i = 0; i < l; i++) {
            if (b[i] > 64 && b[i] < 128) {
                if (cdn[b[i]] > 0) {
                    b[++n] = cdn[b[i]];
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
        cdn[65] = 97;   // A
        cdn[66] = 98;   // B
        cdn[67] = 99;   // C
        cdn[68] = 100;  // D
        cdn[69] = 97;   // dA=E
        cdn[70] = 99;   // dC=F
        cdn[71] = 103;  // G
        cdn[72] = 104;  // H
        cdn[73] = 103;  // I 
        cdn[74] = 103;  // dG=J
        cdn[75] = 107;  // K
        cdn[76] = 116;  // dT=L   
        cdn[77] = 109;  // M
        cdn[78] = 110;  // N
        cdn[82] = 114;  // R
        cdn[83] = 115;  // S
        cdn[84] = 116;  // T
        cdn[85] = 116;  // U
        cdn[86] = 118;  // V
        cdn[87] = 119;  // W
        cdn[89] = 121;  // Y 

        cdn[97] = 97;    // a
        cdn[98] = 98;    // b
        cdn[99] = 99;    // c
        cdn[100] = 100;  // d         
        cdn[101] = 97;   // dA=E
        cdn[102] = 99;   // dC=F
        cdn[103] = 103;  // g
        cdn[104] = 104;  // h  
        cdn[105] = 103;  // i
        cdn[106] = 103;  // dG=J
        cdn[107] = 107;  // k
        cdn[108] = 116;  // dT=L             
        cdn[109] = 109;  // m
        cdn[110] = 110;  // n
        cdn[114] = 114;  // r
        cdn[115] = 115;  // s
        cdn[116] = 116;  // t 
        cdn[117] = 116;  // u
        cdn[118] = 118;  // v
        cdn[119] = 119;  // w
        cdn[121] = 121;  // y  

        int l = source.length();
        byte[] b = source.getBytes();
        int n = -1;
        for (int i = 0; i < l; i++) {
            if (b[i] > 64 && b[i] < 128) {
                if (cdn[b[i]] > 0) {
                    b[++n] = cdn[b[i]];
                }
            }
        }
        return (n > -1) ? (new String(b, 0, n + 1)) : "";
    }
}
