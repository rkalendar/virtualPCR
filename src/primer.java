
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

    public int getLC(String s) {
        return LC(s);
    }

    public int getLC() {
        return LC(seq);
    }

    public int getLen() {
        return lseq;
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

    public double getTm77() {
        return (lseq < 1) ? 0 : (77.1 - dmso + 11.7 * Math.log10(Salt_M) + (41 * (gx + cx + sx + ix + (rx + mx + yx + kx + nx) / 2 + (bx + bx + vx + vx + dx + hx) / 3) - 528) / lseq);
    }

    public double getTm85() {
        return (lseq < 1) ? 0 : (81.5 - dmso + 16.6 * Math.log10(Salt_M) + (41 * (gx + cx + sx + ix + (rx + mx + yx + kx + nx) / 2 + (bx + bx + vx + vx + dx + hx) / 3) - 675) / lseq);
    }

    public double getTm65() {
        return (lseq < 1) ? 0 : (64.9 - dmso + (41 * (gx + cx + sx + ix + (rx + mx + yx + kx + nx) / 2 + (bx + bx + vx + vx + dx + hx) / 3 - 16.4)) / lseq);
    }

    public double getTmS() {
        return (lseq < 1) ? 0 : (2 * (gx + cx + sx + ix + (rx + mx + yx + kx + nx) / 2 + (bx + bx + vx + vx + dx + hx) / 3 + lseq) - dmso);
    }

    public double getMW() {
        double[] d = new double[128];
        d[97] = 313.209;
        d[98] = 307.5293;
        d[99] = 289.184;
        d[100] = 315.5377;
        d[103] = 329.208;
        d[104] = 302.1963;
        d[105] = 314.194;
        d[107] = 316.702;
        d[109] = 301.1965;
        d[110] = 308.9492;
        d[114] = 321.2085;
        d[115] = 309.196;
        d[116] = 304.196;
        d[117] = 290.169;
        d[118] = 310.5337;
        d[119] = 308.7025;
        d[120] = 308.949;
        d[121] = 296.69;
        return (lseq < 1) ? 0 : (ax * d[97] + bx * d[98] + cx * d[99] + dx * d[100] + gx * d[103] + hx * d[104] + ix * d[105] + kx * d[107] + mx * d[109] + nx * d[110] + rx * d[114] + tx * d[116] + sx * d[115] + ux * d[117] + vx * d[118] + wx * d[119] + yx * d[121] - 61.964);
    }

    public double getdH() {
        return dH / 1000d;
    }

    public double getdG() {
        return dG / 1000d;
    }

    public double getdS() {
        return dS;
    }

    public double getGC() {
        return (lseq < 1) ? 0 : ((gx + cx + sx + ix + (rx + mx + yx + kx + nx) / 2 + (bx + bx + vx + vx + dx + hx) / 3) * 100 / lseq);
    }

    public double getA() {
        return (ax + (mx + wx + rx) / 2 + (hx + vx + dx) / 3 + nx / 4);
    }

    public double getT() {
        return (tx + (wx + kx + yx) / 2 + (dx + hx + bx) / 3 + nx / 4);
    }

    public double getG() {
        return (gx + (sx + kx + rx) / 2 + (bx + vx + dx) / 3 + nx / 4);
    }

    public double getC() {
        return (cx + (mx + sx + yx) / 2 + (vx + hx + bx) / 3 + nx / 4);
    }

    public double getI() {
        return ix;
    }

    public double getN() {
        return (bx + dx + hx + kx + mx + nx + rx + sx + vx + wx + yx);
    }

    public double getR() {
        return (ax + ix + (mx + wx + rx + rx + sx + kx + nx) / 2 + (hx + vx + vx + dx + dx + bx) / 3);
    }

    public double getY() {
        return (cx + tx + (mx + sx + yx + wx + kx + yx + nx) / 2 + (vx + hx + bx + dx + hx + bx) / 3);
    }

    public double getDegeneracy() {
        double n = 0;
        if (bx + dx + hx + kx + mx + nx + rx + sx + vx + wx + yx > 0) {
            n = 1;
            if (nx > 0) {
                n = n * nx * 4;
            }
            if (bx > 0) {
                n = n * bx * 3;
            }
            if (dx > 0) {
                n = n * dx * 3;
            }
            if (hx > 0) {
                n = n * hx * 3;
            }
            if (vx > 0) {
                n = n * vx * 3;
            }
            if (kx > 0) {
                n = n * kx * 2;
            }
            if (mx > 0) {
                n = n * mx * 2;
            }
            if (rx > 0) {
                n = n * rx * 2;
            }
            if (sx > 0) {
                n = n * sx * 2;
            }
            if (yx > 0) {
                n = n * yx * 2;
            }
            if (wx > 0) {
                n = n * wx * 2;
            }
        }
        return n;
    }

    public int CountTaq(String s) {
        int k = -0;
        int n = 0;
        int l = s.length();
        if (l == 0) {
            return k;
        }
        for (;;) {
            k = seq.indexOf(s, k);
            if (k == -1) {
                break;
            }
            n++;
            k++;
        }
        return n;
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

        if (seq.equals(dna.ReverseSeq(seq))) {
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
        int r = 0;
        int k = 2;
        int n1 = 0;
        int n2 = 0;
        int n3 = 0;
        int n4 = 0;
        int n5 = 0;
        int n6 = 0;

        int l = source.length();
        if (l < 4) {
            return 100;
        }

        byte[] y1 = new byte[5];
        byte[][] y2 = new byte[5][5];
        byte[][][] y3 = new byte[5][5][5];

        int t = 4 + (l > 16 ? 16 : l - 1);
        if (l > 12) {
            t = t + (l > 65 ? 64 : l - 2);
            k = 3;
        }
        if (l > 48) {
            t = t + (l > 258 ? 256 : l - 3);
            k = 4;
        }
        if (l > 192) {
            t = t + (l > 1027 ? 1024 : l - 4);
            k = 5;
        }
        if (l > 768) {
            t = t + (l > 5000 ? 4096 : l - 5);
            k = 6;
        }

        if (k == 6) {
            byte[][][][] y4 = new byte[5][5][5][5];
            byte[][][][][] y5 = new byte[5][5][5][5][5];
            byte[][][][][][] y6 = new byte[5][5][5][5][5][5];
            for (int i = 0; i < l - 5; i++) {
                n1 = tables.dx[source.charAt(i)];
                n2 = tables.dx[source.charAt(i + 1)];
                n3 = tables.dx[source.charAt(i + 2)];
                n4 = tables.dx[source.charAt(i + 3)];
                n5 = tables.dx[source.charAt(i + 4)];
                n6 = tables.dx[source.charAt(i + 5)];
                if (y1[n1] == 0) {
                    r++;
                    y1[n1] = 1;
                }
                if (y2[n1][n2] == 0) {
                    r++;
                    y2[n1][n2] = 1;
                }
                if (y3[n1][n2][n3] == 0) {
                    r++;
                    y3[n1][n2][n3] = 1;
                }
                if (y4[n1][n2][n3][n4] == 0) {
                    r++;
                    y4[n1][n2][n3][n4] = 1;
                }
                if (y5[n1][n2][n3][n4][n5] == 0) {
                    r++;
                    y5[n1][n2][n3][n4][n5] = 1;
                }
                if (y6[n1][n2][n3][n4][n5][n6] == 0) {
                    r++;
                    y6[n1][n2][n3][n4][n5][n6] = 1;
                }
            }
            if (y5[n2][n3][n4][n5][n6] == 0) {
                r++;
                y5[n2][n3][n4][n5][n6] = 1;
            }
            if (y4[n3][n4][n5][n6] == 0) {
                r++;
                y4[n3][n4][n5][n6] = 1;
            }
            if (y4[n2][n3][n4][n5] == 0) {
                r++;
                y4[n2][n3][n3][n5] = 1;
            }
            if (y3[n4][n5][n6] == 0) {
                r++;
                y3[n4][n5][n6] = 1;
            }
            if (y3[n3][n4][n5] == 0) {
                r++;
                y3[n3][n4][n5] = 1;
            }
            if (y3[n2][n3][n4] == 0) {
                r++;
                y3[n2][n3][n4] = 1;
            }
            if (y2[n5][n6] == 0) {
                r++;
                y2[n5][n6] = 1;
            }
            if (y2[n4][n5] == 0) {
                r++;
                y2[n4][n5] = 1;
            }
            if (y2[n3][n4] == 0) {
                r++;
                y2[n3][n4] = 1;
            }
            if (y2[n2][n3] == 0) {
                r++;
                y2[n2][n3] = 1;
            }
            if (y1[n2] == 0) {
                r++;
                y1[n2] = 1;
            }
            if (y1[n3] == 0) {
                r++;
                y1[n3] = 1;
            }
            if (y1[n4] == 0) {
                r++;
                y1[n4] = 1;
            }
            if (y1[n5] == 0) {
                r++;
                y1[n5] = 1;
            }
            if (y1[n6] == 0) {
                r++;
                y1[n6] = 1;
            }
        }

        if (k == 5) {
            byte[][][][] y4 = new byte[5][5][5][5];
            byte[][][][][] y5 = new byte[5][5][5][5][5];
            for (int i = 0; i < l - 4; i++) {
                n1 = tables.dx[source.charAt(i)];
                n2 = tables.dx[source.charAt(i + 1)];
                n3 = tables.dx[source.charAt(i + 2)];
                n4 = tables.dx[source.charAt(i + 3)];
                n5 = tables.dx[source.charAt(i + 4)];
                if (y1[n1] == 0) {
                    r++;
                    y1[n1] = 1;
                }
                if (y2[n1][n2] == 0) {
                    r++;
                    y2[n1][n2] = 1;
                }
                if (y3[n1][n2][n3] == 0) {
                    r++;
                    y3[n1][n2][n3] = 1;
                }
                if (y4[n1][n2][n3][n4] == 0) {
                    r++;
                    y4[n1][n2][n3][n4] = 1;
                }
                if (y5[n1][n2][n3][n4][n5] == 0) {
                    r++;
                    y5[n1][n2][n3][n4][n5] = 1;
                }
            }
            if (y4[n2][n3][n4][n5] == 0) {
                r++;
                y4[n2][n3][n3][n5] = 1;
            }
            if (y3[n3][n4][n5] == 0) {
                r++;
                y3[n3][n4][n5] = 1;
            }
            if (y3[n2][n3][n4] == 0) {
                r++;
                y3[n2][n3][n4] = 1;
            }
            if (y2[n4][n5] == 0) {
                r++;
                y2[n4][n5] = 1;
            }
            if (y2[n3][n4] == 0) {
                r++;
                y2[n3][n4] = 1;
            }
            if (y2[n2][n3] == 0) {
                r++;
                y2[n2][n3] = 1;
            }
            if (y1[n2] == 0) {
                r++;
                y1[n2] = 1;
            }
            if (y1[n3] == 0) {
                r++;
                y1[n3] = 1;
            }
            if (y1[n4] == 0) {
                r++;
                y1[n4] = 1;
            }
            if (y1[n5] == 0) {
                r++;
                y1[n5] = 1;
            }
        }

        if (k == 4) {
            byte[][][][] y4 = new byte[5][5][5][5];
            for (int i = 0; i < l - 3; i++) {
                n1 = tables.dx[source.charAt(i)];
                n2 = tables.dx[source.charAt(i + 1)];
                n3 = tables.dx[source.charAt(i + 2)];
                n4 = tables.dx[source.charAt(i + 3)];
                if (y1[n1] == 0) {
                    r++;
                    y1[n1] = 1;
                }
                if (y2[n1][n2] == 0) {
                    r++;
                    y2[n1][n2] = 1;
                }
                if (y3[n1][n2][n3] == 0) {
                    r++;
                    y3[n1][n2][n3] = 1;
                }
                if (y4[n1][n2][n3][n4] == 0) {
                    r++;
                    y4[n1][n2][n3][n4] = 1;
                }
            }
            if (y3[n2][n3][n4] == 0) {
                r++;
                y3[n2][n3][n4] = 1;
            }
            if (y2[n3][n4] == 0) {
                r++;
                y2[n3][n4] = 1;
            }
            if (y2[n2][n3] == 0) {
                r++;
                y2[n2][n3] = 1;
            }
            if (y1[n2] == 0) {
                r++;
                y1[n2] = 1;
            }
            if (y1[n3] == 0) {
                r++;
                y1[n3] = 1;
            }
            if (y1[n4] == 0) {
                r++;
                y1[n4] = 1;
            }
        }

        if (k == 3) {
            for (int i = 0; i < l - 2; i++) {
                n1 = tables.dx[source.charAt(i)];
                n2 = tables.dx[source.charAt(i + 1)];
                n3 = tables.dx[source.charAt(i + 2)];
                if (y1[n1] == 0) {
                    r++;
                    y1[n1] = 1;
                }
                if (y2[n1][n2] == 0) {
                    r++;
                    y2[n1][n2] = 1;
                }
                if (y3[n1][n2][n3] == 0) {
                    r++;
                    y3[n1][n2][n3] = 1;
                }
            }
            if (y2[n2][n3] == 0) {
                r++;
                y2[n2][n3] = 1;
            }
            if (y1[n2] == 0) {
                r++;
                y1[n2] = 1;
            }
            if (y1[n3] == 0) {
                r++;
                y1[n3] = 1;
            }
        }

        if (k == 2) {
            for (int i = 0; i < l - 1; i++) {
                n1 = tables.dx[source.charAt(i)];
                n2 = tables.dx[source.charAt(i + 1)];
                if (y1[n1] == 0) {
                    r++;
                    y1[n1] = 1;
                }
                if (y2[n1][n2] == 0) {
                    r++;
                    y2[n1][n2] = 1;
                }
            }
            if (y1[n2] == 0) {
                r++;
                y1[n2] = 1;
            }
        }
        return ((100 * r) / t);
    }
}
