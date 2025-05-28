public final class tools { 

    public static int DNAtest(String source) {
        byte[] cdn = new byte[128];
        cdn[65] = 1;   //A
        cdn[66] = 1;   //B
        cdn[67] = 1;   //C
        cdn[68] = 1;   //D
        cdn[71] = 1;   //G
        cdn[72] = 1;   //H
        cdn[75] = 1;   //K
        cdn[77] = 1;   //M
        cdn[78] = 1;   //N
        cdn[82] = 1;   //R
        cdn[83] = 1;   //S
        cdn[84] = 1;   //T
        cdn[85] = 1;   //U
        cdn[86] = 1;   //V
        cdn[87] = 1;   //W
        cdn[89] = 1;   //Y
        cdn[97] = 1;   //a
        cdn[98] = 1;   //b
        cdn[99] = 1;   //c
        cdn[100] = 1;  //d
        cdn[103] = 1;  //g
        cdn[104] = 1;  //h
        cdn[107] = 1;  //k
        cdn[109] = 1;  //m
        cdn[110] = 1;  //n
        cdn[114] = 1;  //r
        cdn[115] = 1;  //s
        cdn[116] = 1;  //t
        cdn[117] = 1;  //u
        cdn[118] = 1;  //v
        cdn[119] = 1;  //w
        cdn[121] = 1;  //y
        for (int i = 0; i < source.length(); i++) {
            int k = source.charAt(i);
            if (k > 64 && k < 122 && cdn[k] > 0) {
                return 1;
            }
        }
        return 0;
    }

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
        java.text.DecimalFormat f1 = new java.text.DecimalFormat(s);
        return (f1.format(d));
    }

    public static double SeqToDouble(String str) {
        StringBuilder r = new StringBuilder();
        int n = 0;
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
                if (n == 0) {
                    r.append('.');
                    n = 1;
                }
            }
        }
        // return (Double.valueOf("0" + r.toString()).doubleValue());
        return (Double.parseDouble(r.toString()));
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
        return (Integer.parseInt(r.toString()));
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
