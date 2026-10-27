public final class tools { 

    // "Is this a nucleotide character" lookup, previously rebuilt on every call - and rebuilt in
    // full even though the method returns on the first valid character it meets. Read-only.
    // Not tables.cdn: that one maps 'i', this one deliberately does not.
    private static final byte[] ISDNA = new byte[128];

    static {
        ISDNA[65] = 1;   //A
        ISDNA[66] = 1;   //B
        ISDNA[67] = 1;   //C
        ISDNA[68] = 1;   //D
        ISDNA[71] = 1;   //G
        ISDNA[72] = 1;   //H
        ISDNA[75] = 1;   //K
        ISDNA[77] = 1;   //M
        ISDNA[78] = 1;   //N
        ISDNA[82] = 1;   //R
        ISDNA[83] = 1;   //S
        ISDNA[84] = 1;   //T
        ISDNA[85] = 1;   //U
        ISDNA[86] = 1;   //V
        ISDNA[87] = 1;   //W
        ISDNA[89] = 1;   //Y
        ISDNA[97] = 1;   //a
        ISDNA[98] = 1;   //b
        ISDNA[99] = 1;   //c
        ISDNA[100] = 1;  //d
        ISDNA[103] = 1;  //g
        ISDNA[104] = 1;  //h
        ISDNA[107] = 1;  //k
        ISDNA[109] = 1;  //m
        ISDNA[110] = 1;  //n
        ISDNA[114] = 1;  //r
        ISDNA[115] = 1;  //s
        ISDNA[116] = 1;  //t
        ISDNA[117] = 1;  //u
        ISDNA[118] = 1;  //v
        ISDNA[119] = 1;  //w
        ISDNA[121] = 1;  //y
    }

    public static int DNAtest(String source) {
        for (int i = 0; i < source.length(); i++) {
            int k = source.charAt(i);
            if (k > 64 && k < 122 && ISDNA[k] > 0) {
                return 1;
            }
        }
        return 0;
    }


    private static final String DEFAULT_PATTERN = "0.0";
    private static final ThreadLocal<java.text.DecimalFormat> DEFAULT_FORMAT =
            ThreadLocal.withInitial(() -> new java.text.DecimalFormat(DEFAULT_PATTERN,
                    new java.text.DecimalFormatSymbols(java.util.Locale.ROOT)));

    public static String NumberToSeq(double d, int i) {
        String s = "0.0";
        if (i == 0) {
            s = "##";
        }
        if (i > 1) {
            char[] value1 = new char[i];
            for (int j = 0; j < i; j++) {
                value1[j] = '0';
            }
            s = "0." + new String(value1);
        }
        // Force a fixed locale so the decimal separator is always '.', never ',' on
        // comma-locale hosts (ru/kk/de/fr) — the report is tab-delimited and must parse downstream.
        // The one pattern every live call site asks for is cached; DecimalFormat is not
        // thread-safe, hence per-thread rather than one shared instance - a single static would
        // reintroduce exactly the kind of shared mutable state that keeps the engine unparallelisable.
        java.text.DecimalFormat f1 = s.equals(DEFAULT_PATTERN) ? DEFAULT_FORMAT.get()
                : new java.text.DecimalFormat(s, new java.text.DecimalFormatSymbols(java.util.Locale.ROOT));
        return (f1.format(d));
    }

    public static int SeqToInt(String str) {
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

    public static String Strings(int n, char c) {
        if (n < 1) {
            n = 1;
        }
        char[] value1 = new char[n];
        for (int i = 0; i < n; i++) {
            value1[i] = c;
        }
        return new String(value1);
    } 
}
