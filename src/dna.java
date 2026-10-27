public final class dna {

    final public static int maxdna = 1025; // =2^12 or 4^6=4096 (+1); 3^12=531,441; 4^12=16,777,216

    // Complement lookup, shared by AntisenseDNA and ComplementDNA: both built this exact table
    // on every call. Read-only - the loops below mutate the byte[] returned by getBytes(),
    // never the table. Deliberately NOT folded into tables.cdn, which maps different characters.
    private static final byte[] CDNA = new byte[128];
    // Normalisation lookup for DNA(): folds LNA/inosine spellings onto acgt. Also read-only,
    // and also distinct from tables.cdn: i -> g here, and e/f/j/l are mapped here at all.
    private static final byte[] CDN = new byte[128];

    static {
        CDNA[65] = 116;  //A->t
        CDNA[66] = 118;  //B->v
        CDNA[67] = 103;  //C->g
        CDNA[68] = 104;  //D->h
        CDNA[71] = 99;   //G->c
        CDNA[72] = 100;  //H->d
        CDNA[73] = 99;   //I->c
        CDNA[75] = 109;  //K->m
        CDNA[77] = 107;  //M->k
        CDNA[78] = 110;  //N->n
        CDNA[82] = 121;  //R->y
        CDNA[83] = 115;  //S->s
        CDNA[84] = 97;   //T->a
        CDNA[85] = 97;   //U->a
        CDNA[86] = 98;   //V->b
        CDNA[87] = 119;  //W->w
        CDNA[89] = 114;  //Y->r

        CDNA[97] = 116;//  t <- a
        CDNA[98] = 118;//  v <- b
        CDNA[99] = 103;//  g <- c
        CDNA[100] = 104;// h <- d
        CDNA[103] = 99; // c <- g
        CDNA[104] = 100;// d <- h
        CDNA[105] = 99; // i -> c
        CDNA[107] = 109;// m <- k
        CDNA[109] = 107;// k <- m
        CDNA[110] = 110;// n <- n
        CDNA[114] = 121;// y <- r
        CDNA[115] = 115;// s <- s
        CDNA[116] = 97; // a <- t
        CDNA[117] = 97; // a <- u
        CDNA[118] = 98; // b <- v
        CDNA[119] = 119;// w <- w
        CDNA[121] = 114;// r <- y

        CDN[65] = 97;   // A
        CDN[66] = 98;   // B
        CDN[67] = 99;   // C
        CDN[68] = 100;  // D
        CDN[69] = 97;   // dA=E
        CDN[70] = 99;   // dC=F
        CDN[71] = 103;  // G
        CDN[72] = 104;  // H
        CDN[73] = 103;  // I 
        CDN[74] = 103;  // dG=J
        CDN[75] = 107;  // K
        CDN[76] = 116;  // dT=L   
        CDN[77] = 109;  // M
        CDN[78] = 110;  // N
        CDN[82] = 114;  // R
        CDN[83] = 115;  // S
        CDN[84] = 116;  // T
        CDN[85] = 116;  // U
        CDN[86] = 118;  // V
        CDN[87] = 119;  // W
        CDN[89] = 121;  // Y 

        CDN[97] = 97;    // a
        CDN[98] = 98;    // b
        CDN[99] = 99;    // c
        CDN[100] = 100;  // d             
        CDN[101] = 97;   // dA=E
        CDN[102] = 99;   // dC=F
        CDN[103] = 103;  // g
        CDN[104] = 104;  // h  
        CDN[105] = 103;  // i
        CDN[106] = 103;  // dG=J
        CDN[107] = 107;  // k
        CDN[108] = 116;  // dT=L                 
        CDN[109] = 109;  // m
        CDN[110] = 110;  // n
        CDN[114] = 114;  // r
        CDN[115] = 115;  // s
        CDN[116] = 116;  // t 
        CDN[117] = 116;  // u
        CDN[118] = 118;  // v
        CDN[119] = 119;  // w
        CDN[121] = 121;  // y  
    }

    public static String AntisenseDNA(String source) {
        byte[] b = source.getBytes();
        for (int i = 0; i < source.length(); i++) {
            if (CDNA[b[i]] > 0) {
                b[i] = CDNA[b[i]];
            }
        }
        return new String(b);
    }

    public static String ReverseSeq(String source) {
        return new StringBuilder(source).reverse().toString();
    }

    public static String ComplementDNA(String source) {
        byte[] b = source.getBytes();
        int l = source.length();
        int n = l / 2;
        for (int i = 0; i < n; i++) {
            int li = b[i] & 0xFF;
            int ri = b[l - i - 1] & 0xFF;
            // Reverse-complement: complement each end independently (keep unmapped chars verbatim) and always swap.
            byte lc = (li < 128 && CDNA[li] > 0) ? CDNA[li] : b[i];
            byte rc = (ri < 128 && CDNA[ri] > 0) ? CDNA[ri] : b[l - i - 1];
            b[i] = rc;
            b[l - i - 1] = lc;
        }
        if ((l % 2) == 1) {
            int ci = b[n] & 0xFF;
            if (ci < 128 && CDNA[ci] > 0) {
                b[n] = CDNA[ci];
            }
        }
        return new String(b);
    }

    public static String[] ConvertDNA(String source) {
        if (source == null || source.isEmpty()) {
            return null;
        }
        
     //   source="aaaaaaaa";
        
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
 //    r[0] = Arrays.copyOf(f, l);      

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

    public static String DNA(String source) {
        if (source == null || source.isEmpty()) {
            return "";
        }
        int l = source.length();
        byte[] b = source.getBytes();
        int n = -1;
        for (int i = 0; i < l; i++) {
            if (b[i] > 64 && b[i] < 128) {
                if (CDN[b[i]] > 0) {
                    b[++n] = CDN[b[i]];
                }
            }
        }
        return (n > -1) ? (new String(b, 0, n + 1)) : "";
    }
}
