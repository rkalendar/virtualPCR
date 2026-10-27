
public final class primer {
//M=(A/C) R=(A/G) W=(A/T) S=(G/C) Y=(C/T) K=(G/T) V=(A/G/C) H=(A/C/T) D=(A/G/T) B=(C/G/T) N=(A/G/C/T), U=T    
//a   b   c   d   g   h   i   k   m   n   r   s   t   u   v   w   y
//97  98  99  100 103 104 105 107 109 110 114 115 116 117 118 119 121    

    private double mx;
    private double rx;
    private double wx;
    private double sx;
    private double yx;
    private double kx;
    private double vx;
    private double hx;
    private double dx;
    private double bx;
    private double nx;
    private double ux;
    private double ax;
    private double tx;
    private double gx;
    private double cx;
    private double ix;
    private double tm;
    private double dH;
    private double dS;
    private double dG;
    private String seq;
    private int lseq;
    private double dmso = 0;
    private double Mg_M = 0;
    private double Na_M = 0.125d;
    private double Salt_M = 0.125d;
    private double pri = 0.2 / 2000000d;

    public primer(String primer, double K_mM, double mgcl_mM, double dmso0, double primerconc) {
        lseq = primer.length();
        seq = primer;
        dmso = (dmso0 > 0 & dmso0 < 100) ? dmso0 * 0.6d : 0.0;
        Mg_M = mgcl_mM > 0 ? mgcl_mM / 1000 : 0.0;       // mM -> M   
        Na_M = K_mM > 0 ? K_mM / 1000 : 0.001;           // mM -> M 
        Salt_M = K_mM + mgcl_mM > 0.001 ? (120 * Math.sqrt(mgcl_mM) + K_mM) / 1000 : 0.125; // mM calculation
        pri = primerconc > 0.001 ? primerconc / 2000000d : 0.2 / 2000000d;          // 0.25 mkM -> M 

        /* 
4000000 - two different  primers -> if the strands are in equal concentration or by (CA - CB/2) if the strands are at different concentrations, where CA and CB are the concentrations of the more concentrated and less concentrated strands, respectively.
2000000 - primer and genomic DNA
1000000 - self dimer
         */
        if (lseq > 1) {
            Calculation();
        }
    }

    public primer(String primer, int f) {
        lseq = primer.length();
        if (lseq > 1) {
            if (f == 0) {
                seq = primer;
                Calculation();
            } else {
                Calculation(primer);
            }
        }
    }

    public primer(String primer) {
        lseq = primer.length();
        seq = primer;
        if (lseq > 1) {
            Calculation2();
        }
    }

    public primer() {
        Salt_M = 0.11d;
        Mg_M = 0;
        Na_M = 0.125d;
        pri = 0.2d / 2000000d;
    }

    public double getTm(String s) {
        lseq = s.length();
        seq = s;
        if (lseq > 1) {
            Calculation();
        }
        return tm;
    }

    public double getTm() {
        if (tm < 4) {
            return (lseq < 11) ? getTmS() : getTm65();
        }
        return tm;
    }

    public double getTm65() {
        return (lseq < 1) ? 0 : (64.9 - dmso + (41 * (gx + cx + sx + ix + (rx + mx + yx + kx + nx) / 2 + (bx + bx + vx + vx + dx + hx) / 3 - 16.4)) / lseq);
    }

    public double getTmS() {
        return (lseq < 1) ? 0 : (2 * (gx + cx + sx + ix + (rx + mx + yx + kx + nx) / 2 + (bx + bx + vx + vx + dx + hx) / 3 + lseq) - dmso);
    }

    /*
     'a   b   c   d   g   h   i   k   m   n   r   s   t   u   v   w   x   y
     '97  98  99  100 103 104 105 107 109 110 114 115 116 117 118 119 120 121
     'M=(A/C) R=(A/G) W=(A/T) S=(G/C) Y=(C/T) K=(G/T) V=(A/G/C) H=(A/C/T) D=(A/G/T) B=(C/G/T) N=(A/G/C/T), U=T, I            
    
     'SantaLucia J (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. Proceedings of the National Academy of Sciences USA 95 (4):1460-1465
     '     dH   dS    dG
     'AA -7600 -21.3 -1000
     'AC -8400 -22.4 -1440
     'AG -7800 -21.0 -1280
     'AT -7200 -20.4 -880
     'CA -8500 -22.7 -1450
     'CC -8000 -19.9 -1840
     'CG -10600 -27.2 -2170
     'CT -7800 -21.0 -1280
     'GA -8200 -22.2 -1300
     'GC -9800 -24.4 -2240
     'GG -8000 -19.9 -1840
     'GT -8400 -22.4 -1440
     'TA -7200 -21.3 -580
     'TC -8200 -22.2 -1300
     'TG -8500 -22.7 -1450
     'TT -7600 -21.3 -1000
     'IA   2300    4.1
     'IG    100   -2.8
     'init  200   -5.7 1960
     'A/T  2200   +6.9 50
     'sym     0   -1.4 430
     */
    private double[] uTm() {
        double dS1 = 4.1;
        double dH1 = 2300;
        double dG1 = 1030;
        double lg = Math.log(Na_M);
        double[] n = new double[128];
        int n1 = 0;
        int n2 = 0;
        int n3 = 0;

        if (tables.cgv[seq.charAt(0)] == 100) {
            dH1 = 100;
            dS1 = -2.8;
            dG1 = 980;
        }
        if (tables.cgv[seq.charAt(lseq - 1)] == 100) {
            dH1 = dH1 + 100;
            dS1 = dS1 + -2.8;
            dG1 = dG1 + 980;
        } else {
            dH1 = dH1 + 2300;
            dS1 = dS1 + 4.1;
            dG1 = dG1 + 1030;
        }

        if (seq.equals(dna.ComplementDNA(seq))) { // self-complementary, not merely a reversed string
            dS1 = dS1 - 1.4;
            dG1 = dG1 + 400;
        }

        for (int i = 0; i < lseq - 1; i++) {
            n3 = seq.charAt(i);
            n[n3]++;
            n1 = tables.dn[n3];
            n3 = seq.charAt(i + 1);
            n2 = tables.dn[n3];
            dH1 = dH1 + oligoparam.dHUni[n1][n2];
            dS1 = dS1 + oligoparam.dSUni[n1][n2];
            dG1 = dG1 + oligoparam.dGUni[n1][n2];
        }
        n[n3]++;

        double cg = (n[99] + n[103] + n[105] + n[115]) / lseq;
        double t1 = ((1.9872 * Math.log(pri) + dS1) / dH1);
        int r = (int) ((100 * Math.sqrt(Mg_M)) / Na_M);

        if (r < 22) {
            t1 = t1 + (((4.29 * cg - 3.95) / 100000) * lg) + 0.0000094 * lg * lg;
        } else {
            double mg = Math.log(Mg_M);
            double a = 0.0000392;
            double d = 0.0000142;
            double g = 0.0000831;
            if (r < 600) {
                a = a * (0.843 - 0.352 * lg * Math.sqrt(Na_M));
                d = d * (1.279 - 0.00403 * lg - 0.00803 * lg * lg);
                g = g * (0.486 - 0.258 * lg + 0.00525 * lg * lg * lg);
            }
            t1 = t1 + a - 0.00000911 * mg + cg * (0.0000626 + d * mg) + (1 / (2 * (lseq - 1))) * (-0.000482 + 0.000525 * mg + g * mg * mg);
        }
        double tm1 = 1 / t1 - 273.15 - dmso;

        //https://www.google.com/patents/EP1103910A1
        // double tm1 = (dH1 / (1.9872 * Math.log(pri) + dS1)) - 273.15 - dmso;
        dS1 = dS1 + (0.368 * (lseq - 1) * lg);       //dS[Na+] = dS[1M Na+] + 0.368 × N × ln[Na+] where N is equal to the oligonucleotide length minus 1
        dG1 = dG1 - (0.114 * (lseq - 1) * lg);       //dG[Na+] = dG[lM Na+] - 0.114 * N * ln[Na+];
        if (tm1 < 4) {
            tm1 = 4;
        }
        if (tm1 > 95) {
            tm1 = 95;
        }
        n[1] = tm1;
        n[2] = dG1;
        n[3] = dH1;
        n[4] = dS1;
        return n;
    }

    private void Calculation() {
        double r[] = uTm();
        tm = r[1];
        dG = r[2];
        dH = r[3];
        dS = r[4];
        ax = r[97];
        bx = r[98];
        cx = r[99];
        dx = r[100];
        gx = r[103];
        hx = r[104];
        ix = r[105];
        kx = r[107];
        mx = r[109];
        nx = r[110];
        rx = r[114];
        sx = r[115];
        tx = r[116];
        ux = r[117];
        vx = r[118];
        wx = r[119];
        yx = r[121];
    }

    private void Calculation2() {
        double[] r = new double[128];
        for (int i = 0; i < lseq; i++) {
            r[seq.charAt(i)]++;
        }
        ax = r[97];
        bx = r[98];
        cx = r[99];
        dx = r[100];
        gx = r[103];
        hx = r[104];
        ix = r[105];
        kx = r[107];
        mx = r[109];
        nx = r[110];
        rx = r[114];
        sx = r[115];
        tx = r[116];
        ux = r[117];
        vx = r[118];
        wx = r[119];
        yx = r[121];
    }

    private void Calculation(String source) {
        double[] r = new double[128];
        if (source == null || source.isEmpty()) {
            return;
        }
        int[] cdn = new int[128];
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
        int t = 0;
        lseq = 0;
        for (int i = 0; i < l; i++) {
            if (b[i] < 122 && b[i] > 64) {
                t = cdn[b[i]];
                if (t > 0) {
                    lseq++;
                    if (t == 1) {
                        r[b[i]]++;
                    } else {
                        r[b[i] + 32]++;
                    }
                }
            }
        }
        ax = r[97];
        bx = r[98];
        cx = r[99];
        dx = r[100];
        gx = r[103];
        hx = r[104];
        ix = r[105];
        kx = r[107];
        mx = r[109];
        nx = r[110];
        rx = r[114];
        sx = r[115];
        tx = r[116];
        ux = r[117];
        vx = r[118];
        wx = r[119];
        yx = r[121];
    }

       private int LC(String source) {
        byte[] d = new byte[128];
        d[97] = 0;
        d[119] = 0;
        d[101] = 0;
        d[108] = 1;
        d[116] = 1;
        d[117] = 1;
        d[98] = 2;
        d[99] = 2;
        d[102] = 2;
        d[104] = 2;
        d[121] = 2;
        d[109] = 2;
        d[115] = 2;
        d[103] = 3;
        d[105] = 3;
        d[106] = 3;
        d[107] = 3;
        d[114] = 3;
        d[118] = 3;
        d[100] = 3;
        d[110] = 4;

        int l = source.length();
        if (l < 6) {
            return 100;
        }

        int r = 0;
        int y1 = 0;                   // bitmask, 5 bit
        int y2 = 0;                   // bitmask, 25 bit
        long y3a = 0, y3b = 0;        // bitmask, 125 bit (64+61)
        byte[] y4 = new byte[625];    // a flat array instead of [5][5][5][5]

        int t = 4 + Math.min(l - 1, 16) + Math.min(l - 2, 64) + Math.min(l - 3, 256);

        // window initialization
        int n1 = d[source.charAt(0)];
        int n2 = d[source.charAt(1)];
        int n3 = d[source.charAt(2)];

        for (int i = 3; i < l; i++) {
            int n4 = d[source.charAt(i)];

            // y1 — bitmask
            int bit1 = 1 << n1;
            if ((y1 & bit1) == 0) {
                r++;
                y1 |= bit1;
            }

            // y2 — bitmask
            int bit2 = 1 << (n1 * 5 + n2);
            if ((y2 & bit2) == 0) {
                r++;
                y2 |= bit2;
            }

            // y3 — two long masks
            int idx3 = n1 * 25 + n2 * 5 + n3;
            if (idx3 < 64) {
                long m = 1L << idx3;
                if ((y3a & m) == 0) {
                    r++;
                    y3a |= m;
                }
            } else {
                long m = 1L << (idx3 - 64);
                if ((y3b & m) == 0) {
                    r++;
                    y3b |= m;
                }
            }

            // y4 — flat array
            int idx4 = n1 * 125 + n2 * 25 + n3 * 5 + n4;
            if (y4[idx4] == 0) {
                r++;
                y4[idx4] = 1;
            }

            // window shift — one read operation instead of four
            n1 = n2;
            n2 = n3;
            n3 = n4;
        }

        // tail: n1, n2, n3 = the last three characters (positions l-3, l-2, l-1)
        // y3[n1][n2][n3]
        int idx3 = n1 * 25 + n2 * 5 + n3;
        if (idx3 < 64) {
            long m = 1L << idx3;
            if ((y3a & m) == 0) {
                r++;
                y3a |= m;
            }
        } else {
            long m = 1L << (idx3 - 64);
            if ((y3b & m) == 0) {
                r++;
                y3b |= m;
            }
        }

        // y2[n2][n3]
        int bit2 = 1 << (n2 * 5 + n3);
        if ((y2 & bit2) == 0) {
            r++;
            y2 |= bit2;
        }

        // y2[n1][n2]
        bit2 = 1 << (n1 * 5 + n2);
        if ((y2 & bit2) == 0) {
            r++;
            y2 |= bit2;
        }

        // y1[n1], y1[n2], y1[n3]
        int bit1 = 1 << n1;
        if ((y1 & bit1) == 0) {
            r++;
            y1 |= bit1;
        }
        bit1 = 1 << n2;
        if ((y1 & bit1) == 0) {
            r++;
            y1 |= bit1;
        }
        bit1 = 1 << n3;
        if ((y1 & bit1) == 0) {
            r++;
            y1 |= bit1;
        }

        return (100 * r) / t;
    }
}
