
import java.util.HashMap;

public final class InSilicoPCR2 {

    public void SetSequences(String[] seq, String[] sname) {
        this.seq = seq;
        this.sname = sname;
        nseq = seq.length;
    }

    public void SetCurrentFileName(String a) {
        CurrentFileName = a;
    }

    public void SetLookR_Fprimes(boolean a) {
        RFpairs = a;
    }

    public void SetCTbisulfate(boolean a) {
        bis = a;
    }

    public void SetShowPCRProducts(boolean a) {
        pcr_predict = a;
    }

    public void SetShowPrimerAlignment(boolean a) {
        alignment = a;
    }

    public void SetShowPrimerAlignmentPCRproduct(boolean a) {
        PCRmatch_alignment = a;
    }

    public void SetShowOnlyAmplicons(boolean a) {
        ShowOnlyAmplicons = a;
    }

    public void SetShowPCRproductCalculation(boolean a) {
        CalculatePCRproduct = a;
    }

    public void SetProductMaxLength(int v) {
        mxPCRsize = (v < 50000) ? v - 1 : 50000;
    }

    public void SetProductMinLength(int v) {
        mnPCRsize = (v < 12) ? 12 : v - 1;
    }

    public StringBuilder getResult() {
        return sr;
    }

    public void SetPrimers(String[] primerslist, String[] primersname, String[] oriprimers, boolean is_probe, boolean is_pattern, boolean is_circle, boolean Extract, int err3) {
        probe = is_probe;
        patsearch = is_pattern;
        primers = primerslist;
        pname = primersname;
        nprimers = primerslist.length;
        oprimers = oriprimers;
        SequenceExtract = Extract;
        hashkmer = 100;
        kmer = 12;
        er3 = err3;
        if (er3 < 0) {
            er3 = 0;
        }
        if (er3 > 10) {
            er3 = 10;
        }
        if (er3 > 0) {
            kmer = 9;
        }
        ct = is_circle ? kmer : 0;
        if (!patsearch) {
            LoadPrimers1();
            if (probe) {
                PreparedPrimers(-1);
            } else {
                PreparedPrimers(4);
            }
        } else {
            LoadPrimers2(oriprimers);
            PreparedPrimers(-1);
        }
    }

    public void Run() {
        sr = new StringBuilder();
        FastStart(er3);
    }

    private void LoadPrimers1() {
        pr = new String[nprimers];
        apr = new String[nprimers];
        rpr = new String[nprimers];
        cpr = new String[nprimers];
        pa = new int[nprimers];
        prl = new int[nprimers];
        prtm = new int[nprimers];
        prnk = new int[nprimers];
        tprl = 0;

        for (int i = 1; i < nprimers; i++) {
            primer prm = new primer(primers[i], 50, 0, 0, 0.25d);
            prtm[i] = (int) prm.getTm();
            int l = primers[i].length();
            pr[i] = primers[i];
            if (l < kmer) {
                if (l > kmer - 3) {
                    pa[i] = kmer - l;
                    pr[i] = primers[i];
                    l = kmer;
                    pr[i] = tools.Strings(pa[i], 'n') + primers[i];
                }
            }
            rpr[i] = dna.ReverseSeq(pr[i]);
            apr[i] = dna.AntisenseDNA(pr[i]);
            cpr[i] = dna.ReverseSeq(apr[i]);
            prl[i] = l;
            tprl = tprl + l;
        }

        if (RFpairs) {
            int k = 0;
            for (int i = 1; i < nprimers - 1; i++) {
                if (prnk[i] == 0) {
                    prnk[i] = ++k;
                    for (int j = i + 1; j < nprimers; j++) {
                        if (pname[i].trim().equals(pname[j].trim())) {
                            prnk[j] = prnk[i];
                        }
                    }
                }
            }
        }
    }

    private void LoadPrimers2(String[] prms) {
        pr = new String[nprimers];
        apr = new String[nprimers];
        rpr = new String[nprimers];
        cpr = new String[nprimers];
        pa = new int[nprimers];
        prl = new int[nprimers];
        prtm = new int[nprimers];
        prnk = new int[nprimers];
        tprl = 0;
        int h = 0;
        int z = 0;
        int k = -1;
        int n = 0;
        int l = 0;

        prn = new int[nprimers]; //prn = skolko vsego elements in pattern       
        prx = new int[nprimers]; //ot kakogo primera, proishodi etot primer   
        prxn = new int[nprimers];//n elements in current pattern 1....n       
        ptx1 = new int[nprimers];//pattern location betwenn x1       
        ptx2 = new int[nprimers];//pattern location to x2    

        for (int j = 1; j < nprimers; j++) {
            String[] s = prms[j].split("\\[");//http://docs.oracle.com/javase/6/docs/api/java/util/regex/Pattern.html#sum  "(?<=\\p{Punct})"
            h = s.length;
            z = 0;
            if (h > 0) {
                for (int i = 0; i < h; i++) {
                    if (!isNumeric(s[i])) {
                        if (tools.DNAtest(s[i]) == 1) {
                            if (s[i].contains("@")) {
                                s[i] = dna.ComplementDNA(dna.DNA(s[i]));
                            } else {
                                s[i] = dna.DNA(s[i]);
                            }

                            l = s[i].length();
                            if (l > 4) {
                                if (l < kmer) {
                                    kmer = l;
                                }
                            }
                        }
                    }
                }
            } else {
                if (!isNumeric(s[0])) {
                    if (tools.DNAtest(s[0]) == 1) {
                        if (s[0].contains("@")) {
                            s[0] = dna.ComplementDNA(dna.DNA(s[0]));
                        } else {
                            s[0] = dna.DNA(s[0]);
                        }
                        l = s[0].length();
                        if (l > 4) {
                            if (l < kmer) {
                                kmer = l;
                            }
                        }
                    }
                }
            }
        }

//gccRggaRR[1-100]tagatgggagaaaggtttcccctct[1-50]gaggggagagagcccctccgagattggatc
        for (int j = 1; j < nprimers; j++) {
            String[] s = prms[j].split("\\p{Punct}");//http://docs.oracle.com/javase/6/docs/api/java/util/regex/Pattern.html#sum
            //      String[] s = prms[j].split("(?<=\\p{Punct})");//http://docs.oracle.com/javase/6/docs/api/java/util/regex/Pattern.html#sum  "(?<=\\p{Punct})"
            h = s.length;
            z = 0;
            n = 0;
            if (h > 0) {
                for (int i = 0; i < h; i++) {
                    if (isNumeric(s[i])) {
                        if (z == 2) {
                            z = tools.SeqToInt(s[i]);
                            if (z > ptx2[k]) {
                                ptx1[k] = ptx2[k];
                            }
                            ptx2[k] = z;
                            z = 0;
                        }
                        if (z == 1) {
                            ptx2[k] = tools.SeqToInt(s[i]);
                            ptx1[k] = 1;
                            z = 2;
                        }
                    } else {
                        if (tools.DNAtest(s[i]) == 1) {
                            s[i] = dna.DNA(s[i]);
                            l = s[i].length();
                            if (l > 4) {
                                k++;
                                n++;
                                z = 1;
                                prn[j]++;   // current primers containg n elements
                                FillPrimers(k, l, n, s[i]);
                                prx[k] = j; //ot kakogo primera in pattern, proishodi etot primer      
                            }
                        }
                    }
                }
            } else// if (h == 0) 
            {
                prms[j] = dna.DNA(prms[j]);
                l = prms[j].length();
                if (l > 4) {
                    k++;
                    n++;
                    prn[j]++;
                    FillPrimers(k, l, n, prms[j]);
                    prx[k] = j;
                }
            }
        }

        nprimers = k + 1;
        if (RFpairs) {
            k = 0;
            for (int i = 1; i < nprimers - 1; i++) {
                int l1 = pname[i].trim().length() - 1;
                if (l1 > 0 && prnk[i] == 0) {
                    k++;
                    prnk[i] = k;
                    String n1 = pname[i].trim().substring(0, l1);
                    for (int j = i + 1; j < nprimers; j++) {
                        int l2 = pname[j].trim().length() - 1;
                        if (l2 > 0 && prnk[j] == 0) {
                            if (pname[i].trim().substring(0, l2).equals(n1)) {
                                prnk[j] = prnk[i];
                            }
                        }
                    }
                }
            }
        }

    }

    private void FillPrimers(int k, int l, int n, String s) {
        String[] a;
        int[] t;

        if (k >= pr.length) {
            a = new String[k + nprimers];
            System.arraycopy(pr, 0, a, 0, k);
            pr = a;

            a = new String[k + nprimers];
            System.arraycopy(apr, 0, a, 0, k);
            apr = a;

            a = new String[k + nprimers];
            System.arraycopy(rpr, 0, a, 0, k);
            rpr = a;

            a = new String[k + nprimers];
            System.arraycopy(cpr, 0, a, 0, k);
            cpr = a;

            t = new int[k + nprimers];
            System.arraycopy(prxn, 0, t, 0, k);
            prxn = t;

            t = new int[k + nprimers];
            System.arraycopy(prx, 0, t, 0, k);
            prx = t;

            t = new int[k + nprimers];
            System.arraycopy(ptx1, 0, t, 0, k);
            ptx1 = t;

            t = new int[k + nprimers];
            System.arraycopy(ptx2, 0, t, 0, k);
            ptx2 = t;

            t = new int[k + nprimers];
            System.arraycopy(pa, 0, t, 0, k);
            pa = t;

            t = new int[k + nprimers];
            System.arraycopy(prl, 0, t, 0, k);
            prl = t;

            t = new int[k + nprimers];
            System.arraycopy(prtm, 0, t, 0, k);
            prtm = t;

            t = new int[k + nprimers];
            System.arraycopy(prnk, 0, t, 0, k);
            prnk = t;
        }
        primer prm = new primer(s, 0.05d, 0, 0, 0.25d);
        prtm[k] = (int) prm.getTm();
        prxn[k] = n;
        pr[k] = s;
        if (l < kmer) {
            if (l > kmer - 3) {
                pa[k] = kmer - l;
                pr[k] = s;
                l = kmer;
                pr[k] = tools.Strings(pa[k], 'n') + s;
            }
        }
        rpr[k] = dna.ReverseSeq(pr[k]);
        apr[k] = dna.AntisenseDNA(pr[k]);
        cpr[k] = dna.ReverseSeq(apr[k]);
        prl[k] = l;
        tprl = tprl + l;
    }

    private void PreparedPrimers(int z) {
        int i = 0;
        int x = 0;
        int b = 0;
        int p = 0;
        int d = 0;
        int t = 0;
        int w = -1;

        int lp = 0;
        int n0 = 0;
        int n1 = 0;
        int na = 0;

        String[] R1;
        String s = "";

        px = new HashMap<>();
        k1 = new int[tprl + tprl + tprl + tprl]; //n elements 

        for (i = 1; i < nprimers; i++) {
            lp = pr[i].length();
            if (lp > kmer - 1) {
                d = z;
                if (d < 0 || d > lp - kmer) {
                    d = 1 + lp - kmer;
                }

                for (b = 0; b < d; b++) {
                    x = lp - kmer - b;
                    na = 1;
                    for (t = 0; t < kmer; t++) {
                        na = na * tables.t2[pr[i].charAt(x + t)][0];
                    }
                    if (na < dna.maxdna) {
                        if (na == 1) {
                            n0++;
                            s = pr[i].substring(x, x + kmer);
                            if (px.containsKey(s)) {
                                p = px.get(s);
                                k1[p]++;
                            } else {
                                w++;
                                px.put(s, w);
                                if (w >= k1.length) {
                                    int[] kc = new int[w + tprl];
                                    System.arraycopy(k1, 0, kc, 0, w);
                                    k1 = kc;
                                }
                                k1[w]++;
                            }
                        } else {
                            R1 = dna.ConvertDNA(pr[i].substring(x, x + kmer));
                            for (String R11 : R1) {
                                n0++;
                                if (px.containsKey(R11)) {
                                    p = px.get(R11);
                                    k1[p]++;
                                } else {
                                    w++;
                                    px.put(R11, w);
                                    if (w >= k1.length) {
                                        int[] kc = new int[w + tprl];
                                        System.arraycopy(k1, 0, kc, 0, w);
                                        k1 = kc;
                                    }
                                    k1[w]++;
                                }
                            }
                        }
                    }
// reverse direction                    
                    x = b;
                    na = 1;
                    for (t = 0; t < kmer; t++) {
                        na = na * tables.t2[cpr[i].charAt(x + t)][0];
                    }
                    if (na < dna.maxdna) {
                        if (na == 1) {
                            n0++;
                            s = cpr[i].substring(x, x + kmer);
                            if (px.containsKey(s)) {
                                p = px.get(s);
                                k1[p]++;
                            } else {
                                w++;
                                px.put(s, w);
                                if (w >= k1.length) {
                                    int[] kc = new int[w + tprl];
                                    System.arraycopy(k1, 0, kc, 0, w);
                                    k1 = kc;
                                }
                                k1[w]++;
                            }
                        } else {
                            R1 = dna.ConvertDNA(cpr[i].substring(x, x + kmer));
                            for (String R11 : R1) {
                                n0++;
                                if (px.containsKey(R11)) {
                                    p = px.get(R11);
                                    k1[p]++;
                                } else {
                                    w++;
                                    px.put(R11, w);
                                    if (w >= k1.length) {
                                        int[] kc = new int[w + tprl];
                                        System.arraycopy(k1, 0, kc, 0, w);
                                        k1 = kc;
                                    }
                                    k1[w]++;
                                }
                            }
                        }
                    }

                }
            }
        }

        u1 = new int[n0];//locations at ->k3       
        u2 = new int[n0];//real locations x, primer            
        p1 = new int[n0];//primer index         
        w = 0;
        for (i = 1; i < nprimers; i++) {
            lp = pr[i].length();
            if (lp > kmer - 1) {
                d = z;
                if (d < 0 || d > lp - kmer) {
                    d = 1 + lp - kmer;
                }
                for (b = 0; b < d; b++) {
                    x = lp - kmer - b;
                    na = 1;
                    for (t = 0; t < kmer; t++) {
                        na = na * tables.t2[pr[i].charAt(x + t)][0];
                    }
                    if (na < dna.maxdna) {
                        if (na == 1) {
                            s = pr[i].substring(x, x + kmer);
                            p = px.get(s);
                            if (k1[p] > 0) {
                                u1[p] = w;
                                p1[w] = i;
                                u2[w] = x;
                                w = w + k1[p];
                                k1[p] = -1;
                            } else {
                                n1 = u1[p] - k1[p];
                                u2[n1] = x;
                                p1[n1] = i;
                                k1[p]--; // n -elements
                            }
                        } else {
                            R1 = dna.ConvertDNA(pr[i].substring(x, x + kmer));
                            for (String R11 : R1) {
                                p = px.get(R11);
                                if (k1[p] > 0) {
                                    u1[p] = w;
                                    p1[w] = i;
                                    u2[w] = x;
                                    w = w + k1[p];
                                    k1[p] = -1;
                                } else {
                                    n1 = u1[p] - k1[p];
                                    u2[n1] = x;
                                    p1[n1] = i;
                                    k1[p]--; // n -elements
                                }
                            }
                        }
                    }

// antisense direction                    
                    x = b;
                    na = 1;
                    for (t = 0; t < kmer; t++) {
                        na = na * tables.t2[cpr[i].charAt(x + t)][0];
                    }
                    if (na < dna.maxdna) {
                        if (na == 1) {
                            s = cpr[i].substring(x, x + kmer);
                            p = px.get(s);
                            if (k1[p] > 0) {
                                u1[p] = w;
                                u2[w] = -(x + 1);//0 -> -1 problem 0
                                p1[w] = i;
                                w = w + k1[p];
                                k1[p] = -1;
                            } else {
                                n1 = u1[p] - k1[p];
                                u2[n1] = -(x + 1);//0 -> -1 problem 0
                                p1[n1] = i;
                                k1[p]--; // n -elements
                            }
                        } else {
                            R1 = dna.ConvertDNA(cpr[i].substring(x, x + kmer));
                            for (String R11 : R1) {
                                p = px.get(R11);
                                if (k1[p] > 0) {
                                    u1[p] = w;
                                    p1[w] = i;
                                    u2[w] = -(x + 1);//0 -> -1 problem 0
                                    w = w + k1[p];
                                    k1[p] = -1;
                                } else {
                                    n1 = u1[p] - k1[p];
                                    u2[n1] = -(x + 1);//0 -> -1 problem 0
                                    p1[n1] = i;
                                    k1[p]--; // n -elements
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private void AddRecord(int x, int j, int p, int l1, int Err3, LinkList xr) {
        char c1;
        char c2;
        int r = 0;
        int l = 0;
        int t = 0;
        int e = 0;
        int x3 = 0;
        int x5 = 0;

        for (int g = 0; g < -k1[p]; g++) {
            int y = u2[u1[p] + g];
            int d = p1[u1[p] + g]; // primer index
            int l2 = pr[d].length();
            int v = l2 < 22 ? l2 - kmer - 2 : 20 - kmer; // test
            v = v - Err3;

            if (y > -1) // direct Forward
            {
                x3 = x - y + l2 - 1;
                x5 = x3 - l2 + 1;
                r = 0;
                l = 0;
                if (x3 > ct - 1) {

                    if (pvx1[d] < x3) {
                        e = 1;
                        for (t = kmer; t < l2 - y; t++) {
                            if (ct > 0) {
                                if (x + t + 1 > l1)// cicle
                                {
                                    c2 = seq[j].charAt(x + t - l1);
                                } else {
                                    c2 = seq[j].charAt(x + t);
                                }
                            } else {
                                if (x + t + 1 > l1) {
                                    break;
                                }
                                c2 = seq[j].charAt(x + t);
                            }
                            c1 = pr[d].charAt(y + t);
                            if (tables.aDNA[c1][tables.cdna[c2]] < tables.smth) { //<= tables.smth
                                r++;
                                if (r > er3 || e == 0) {
                                    break;
                                }
                                e = 0;
                            } else {
                                e = 1;
                            }
                            l++;
                        }
                        if (r <= Err3) {
                            e = 1;
                            for (t = 1; x - t > 0 && y - t > 0; t++) {
                                c1 = pr[d].charAt(y - t);
                                c2 = seq[j].charAt(x - t);
                                if (tables.aDNA[c1][tables.cdna[c2]] < tables.smth) { //<= tables.smmax
                                    r++;
                                    if (r > er3 || e == 0) {
                                        break;
                                    }
                                    e = 0;
                                } else {
                                    e = 1;
                                }
                                l++;
                            }
                        }
                        if (l > v & l > 0) {
                            pvx1[d] = x3;
                            xr.insert(x3, x5, d, 0);
                        }
                    }
                } //ct
            } else // reverse
            {
                y = -y - 1;// -1->0 solving problem with 0
                x3 = x - y;
                x5 = x - y + l2 - 1;
                r = 0;
                l = 0;
                if (x3 + Err3 + 1 > ct) {

                    if (pvx2[d] < x3) {
                        e = 1;
                        for (t = 1; t < y + 1; t++) {
                            if (x - t < 0) {
                                break;
                            }
                            c1 = rpr[d].charAt(y - t);
                            c2 = seq[j].charAt(x - t);
                            if (tables.aDNA[c1][c2] < tables.smth) {
                                r++;
                                if (r > er3 || e == 0) {
                                    break;
                                }
                                e = 0;
                            } else {
                                e = 1;
                            }
                            l++;
                        }
                        if (r <= Err3) {
                            e = 1;
                            for (t = kmer; y + t < l2; t++) {
                                if (ct > 0) {
                                    if (x + t > l1 - 1)// cicle
                                    {
                                        c2 = seq[j].charAt(x + t - l1);
                                    } else {
                                        c2 = seq[j].charAt(x + t);
                                    }
                                } else {
                                    if (x + t > l1 - 1) {
                                        break;
                                    }
                                    c2 = seq[j].charAt(x + t);
                                }
                                c1 = rpr[d].charAt(y + t);
                                if (tables.aDNA[c1][c2] < tables.smth) {
                                    r++;
                                    if (r > er3 || e == 0) {
                                        break;
                                    }
                                    e = 0;
                                } else {
                                    e = 1;
                                }
                                l++;
                            }
                        }
                        if (l > v & l > 0) {
                            pvx2[d] = x3;
                            xr.insert(x3, x5, d, 1);
                        }
                    }
                } //ct                                    
            }
        }
    }

    private void FastStart(int Err3) {
        int x = 0;
        int t = 0;
        int p = 0;
        int na = 0;
        int l1 = 0;
        byte[] ax = new byte[kmer];

        String[] R1;
        String s = "";
        for (int j = 0; j < nseq; j++) {
            nsites = 0;
            l1 = seq[j].length();
            if (l1 > kmer) {

                LinkList xr = new LinkList();

                pvx1 = new int[nprimers];
                pvx2 = new int[nprimers];
                for (t = 0; t < nprimers; t++) {
                    pvx1[t] = -1;
                    pvx2[t] = -1;
                }
                na = 1;
                for (t = 0; t < kmer - 1; t++) {
                    ax[t] = tables.t3[seq[j].charAt(t)][0];
                    na = na * ax[t];
                }
                for (int i = kmer - 1; i < l1 + ct; i++) {
                    x = 1 + i - kmer;
                    if (i > l1 - 1) { // circle, ct=kmer
                        ax[kmer - 1] = tables.t3[seq[j].charAt(i - l1)][0];
                        if (x < l1) {
                            s = seq[j].substring(x) + seq[j].substring(0, i - l1 + 1);
                        } else {
                            s = seq[j].substring(i - l1, i - l1 + kmer);
                        }
                    } else {
                        ax[kmer - 1] = tables.t3[seq[j].charAt(i)][0];
                        s = seq[j].substring(x, x + kmer);
                    }
                    na = na * ax[kmer - 1];
                    if (na < dna.maxdna) {
                        if (na == 1) {
                            if (!px.containsKey(s)) {
                                continue;
                            }
                            p = px.get(s);
                            AddRecord(x, j, p, l1, Err3, xr);
                        } else {
                            if ((ax[kmer - 1] * ax[kmer - 2] * ax[kmer - 3]) < 64 && (ax[0] * ax[1] * ax[2]) < 64) {
                                R1 = dna.ConvertDNA(seq[j].substring(x, x + kmer));
                                for (String R11 : R1) {
                                    if (!px.containsKey(R11)) {
                                        continue;
                                    }
                                    p = px.get(R11);
                                    AddRecord(x, j, p, l1, Err3, xr);
                                }
                            }
                        }
                    }
                    na = na / ax[0];
                    for (t = 0; t < kmer - 1; t++) {
                        ax[t] = ax[t + 1];
                    }
                }
                if (!xr.isEmpty()) // result
                {
                    if (patsearch) {
                        xr.shellSort();
                        PrintResult2(xr, j);
                        sr.append("\nAmount of sites: ").append(nsites).append("\n\n");

                    } else {
                        xr.shellSort();
                        PrintResult1(xr, j);
                    }
                }

            }
        }
    }

    private String Xconvertor(int x5, int x3, int l1, int ct, int drt) {
        int n2 = x5 + 1;
        if (x5 < 1) {
            n2 = 1;
            if (ct > 0) {
                n2 = 1 + l1 + x5;
            }
        }
        if (x5 > l1 - 1) {
            n2 = l1;
            if (ct > 0) {
                n2 = 1 + x5 - l1;
            }
        }
        int n3 = x3 + 1;
        if (x3 < 1) {
            n3 = 1;
            if (ct > 0) {
                n3 = 1 + l1 + x3;
            }
        }
        if (x3 > l1 - 1) {
            n3 = l1;
            if (ct > 0) {
                n3 = 1 + x3 - l1;
            }
        }

        return switch (drt) {
            case 0 ->
                n2 + "->" + n3;
            case 1 ->
                n3 + "<-" + n2;
            default ->
                n2 + "<-" + n3;
        };

    }

    private void PrintPrimerAlignment(int op, int n, int l, int x5, int x3, int pn) {
        int y = 0;
        int x = 0;
        int t = 0;
        int l1 = seq[n].length();
        int n1 = 0;
        int n3 = 0;
        double tm = 0;

        String s = "";
        String s1 = "";
        String s2 = "";
        char[] value;

        if (op == 0) {
            x = x5;
            y = x3;
            s1 = "";
            s2 = "";
            if (x < 0) {
                s1 = (ct > 0) ? seq[n].substring(-x) : tools.Strings(-x, ' ');
                x = 0;
            }

            if (y > l1 - 1) {
                s2 = (ct > 0) ? seq[n].substring(0, 1 + y - l1) : tools.Strings(1 + y - l1, ' ');
                y = l1 - 1;
            }
            s = s1 + dna.AntisenseDNA(seq[n].substring(x, y + 1)) + s2;
            tm = Melting.Tmelting(pr[pn], s, 0.2d, 55d, 1, 0, 4);

            value = new char[l + 2];
            value[0] = ' ';
            value[1] = ' ';

            n1 = 0;
            for (t = 0; t < l; t++) {
                value[t + 2] = ' ';
                n3 = tables.aDNA[pr[pn].charAt(t)][s.charAt(t)];
                if (n3 > tables.smmin) {
                    value[t + 2] = n3 > tables.smth ? '|' : ':';
                }
                n1 = n1 + n3;
            }
            n1 = n1 / l;

            x = x5 - 2;
            y = x3 + 4;
            s1 = "";
            s2 = "";
            if (x < 0) {
                s1 = (ct > 0) ? seq[n].substring(-x) : tools.Strings(-x, ' ');
                x = 0;
            }
            if (y > l1 - 1) {
                s2 = (ct > 0) ? seq[n].substring(0, 1 + y - l1) : tools.Strings(1 + y - l1, ' ');
                y = l1 - 1;
            }
            s = dna.AntisenseDNA(s1 + seq[n].substring(x, y + 1) + s2);
            s1 = (pa[pn] > 0) ? tools.Strings(pa[pn], ' ') + "5-" : "5-";
            sr.append("Position: ").append(Xconvertor(x5, x3, l1, ct, 0)).append("\t").append(n1).append("%\tTm=").append(tools.NumberToSeq(tm, 1)).append("°C\n");
            sr.append(s1).append(primers[pn]).append("->\n").append(new String(value)).append("\n").append(s).append("\n");
        } else ///////////////// reverse primer
        {
            x = x3;
            y = x5;
            s1 = "";
            s2 = "";
            if (x < 0) {
                s1 = ct > 0 ? seq[n].substring(l1 + x) : tools.Strings(-x, ' ');
                x = 0;
            }

            if (y > l1 - 1) {
                s2 = ct > 0 ? seq[n].substring(0, 1 + y - l1) : tools.Strings(1 + y - l1, ' ');
                y = l1 - 1;
            }
            s = s1 + seq[n].substring(x, y + 1) + s2;
            tm = Melting.Tmelting(rpr[pn], s, 0.2d, 55d, 1, 0, 4);

            value = new char[l + 2];
            value[0] = ' ';
            value[1] = ' ';

            n1 = 0;
            for (t = 0; t < l; t++) {
                value[t + 2] = ' ';
                n3 = tables.aDNA[rpr[pn].charAt(t)][s.charAt(t)];
                if (n3 > tables.smmin) {
                    value[t + 2] = n3 > tables.smth ? '|' : ':';
                }
                n1 = n1 + n3;
            }
            n1 = n1 / l;

            x = x3 - 2;
            y = x5 + 4;
            s1 = "";
            s2 = "";
            if (x < 0) {
                s1 = (ct > 0) ? seq[n].substring(l1 + x) : tools.Strings(-x, ' ');
                x = 0;
            }
            if (y > l1 - 1) {
                s2 = (ct > 0) ? seq[n].substring(0, 1 + y - l1) : tools.Strings(1 + y - l1, ' ');
                y = l1 - 1;
            }
            s = s1 + seq[n].substring(x, y + 1) + s2;
            sr.append("Position: ").append(Xconvertor(x5, x3, l1, ct, 1)).append("\t").append(n1).append("%\tTm=").append(tools.NumberToSeq(tm, 1)).append("°C\n");
            sr.append("<-").append(rpr[pn]).append("-5\n").append(new String(value)).append("\n").append(s).append("\n");
        }
    }

    private void PrintResult1(LinkList lnk, int n) {
        int i = 0;
        int j = 0;
        int x0 = 0;
        int x1 = 0;
        int x2 = 0;
        int nampls = 0;

        int k = lnk.Amount();
        int l1 = seq[n].length();
        int n0 = -1;

        String s = "";
        String s1 = "";
        String s2 = "";

        int[] op = lnk.getOriantation();
        int[] x3 = lnk.getX3();
        int[] x5 = lnk.getX5();
        int[] pn = lnk.getPrimers();

        java.util.Map<Integer, Integer> map = new java.util.TreeMap<>(java.util.Collections.reverseOrder());
        java.util.Map<Integer, Integer> map2 = new java.util.TreeMap<>();
        if (!lnk.isEmpty()) {
            sr.append("__________________________________________________\n In silico Primer(s) search for: ").append(CurrentFileName).append("//").append(sname[n]).append(":\n\n");
        }
        if (alignment) {
            for (i = 0; i < k; i++) {
                if (n0 != pn[i]) {
                    n0 = pn[i];
                    sr.append(pname[pn[i]]).append("\n").append(primers[pn[i]]).append("\n\n");
                }
                PrintPrimerAlignment(op[i], n, prl[pn[i]], x5[i], x3[i], pn[i]);
                sr.append("\n");
            }
        }

        if (pcr_predict) {
            for (i = 0; i < k - 1; i++) {
                for (j = i + 1; j < k; j++) {
                    if ((op[i] != op[j]) && prnk[pn[i]] == prnk[pn[j]]) {
                        int PCRproduct1 = 0;
                        int PCRproduct2 = 0;
                        x0 = 0;
                        if (op[i] == 0) { //  Fi-> <-jR
                            x0 = 1;
                            x1 = x5[i];
                            x2 = x5[j];
                            PCRproduct1 = 1 + x5[j] - x5[i]; //<-jR Fi->  F(i)->L1___1<-(j) ' Fi-> <-jR   
                            if (ct > 0) {
                                if (x5[i] > x5[j]) {
                                    PCRproduct2 = l1 - x5[i] + x5[j];
                                } else if (x5[j] > x5[i]) {
                                    PCRproduct2 = l1 + x5[j] - x5[i];
                                }
                            }
                        } else {  // Fj-> <-iR
                            x0 = 2;
                            x1 = x5[j];
                            x2 = x5[i];
                            PCRproduct1 = 1 + x5[i] - x5[j]; //<-iR Fj->  F(j)->L1---1<-(i) ' Fj-> <-iR                    
                            if (ct > 0) {
                                if (x5[j] > x5[i]) {
                                    PCRproduct2 = l1 - x5[j] + x5[i];
                                } else if (x5[i] > x5[j]) {
                                    PCRproduct2 = l1 + x5[i] - x5[j];
                                }
                            }
                        }

                        if ((PCRproduct1 > mnPCRsize && PCRproduct1 < mxPCRsize) | (PCRproduct2 > mnPCRsize && PCRproduct2 < mxPCRsize)) {
                            if (x0 == 1) {
                                s1 = Xconvertor(x5[i], x3[i], l1, ct, 0);
                                s2 = Xconvertor(x5[j], x3[j], l1, ct, 1);
                            }
                            if (x0 == 2) {
                                s1 = Xconvertor(x5[i], x3[i], l1, ct, 1);
                                s2 = Xconvertor(x5[j], x3[j], l1, ct, 0);
                            }
                            if (x1 < 0) {
                                x1 = 0;
                            }
                            if (x1 > l1) {
                                x1 = l1;
                            }
                            if (x2 < 0) {
                                x2 = 0;
                            }
                            if (x2 > l1) {
                                x2 = l1;
                            }

                            if (SequenceExtract) {
                                if (PCRproduct1 > mnPCRsize && PCRproduct1 < mxPCRsize) {
                                    boolean f1 = false;
                                    if (map2.containsKey(x1)) {
                                        if (map2.get(x1) != x2) {
                                            f1 = true;
                                        }
                                    } else {
                                        f1 = true;
                                    }
                                    if (f1) {
                                        map2.put(x1, x2);
                                        sr.append("\n>").append(x1 + 1).append("-").append(x2 + 1).append(" Amplicon size: ").append(PCRproduct1).append("bp\t").append("Ta=").append(tools.NumberToSeq(Melting.Ta(PCRproduct1, prtm[pn[i]], prtm[pn[j]]), 0)).append("°C\n");
                                        sr.append(seq[n].substring(x1, x2)).append("\n");
                                    }

                                    if (PCRproduct2 > mnPCRsize && PCRproduct2 < mxPCRsize) {
                                        if (x0 == 1) {
                                            sr.append("\n>").append(x1 + 1).append("-").append(x2 + 1).append(" Amplicon size: ").append(PCRproduct2).append("bp\t").append("Ta=").append(tools.NumberToSeq(Melting.Ta(PCRproduct2, prtm[pn[i]], prtm[pn[j]]), 0)).append("°C\n");
                                            sr.append(seq[n].substring(x1)).append(seq[n].substring(0, x2)).append("\n");

                                        } else if (x0 == 2) {
                                            sr.append("\n>").append(x1 + 1).append("-").append(x2 + 1).append(" Amplicon size: ").append(PCRproduct2).append("bp\t").append("Ta=").append(tools.NumberToSeq(Melting.Ta(PCRproduct2, prtm[pn[i]], prtm[pn[j]]), 0)).append("°C\n");
                                            sr.append(seq[n].substring(x2)).append(seq[n].substring(0, x1)).append("\n");
                                        }
                                    }
                                }
                            }

                            if (ShowOnlyAmplicons) {
                                if (PCRproduct1 > mnPCRsize && PCRproduct1 < mxPCRsize) {
                                    map.put(PCRproduct1, Melting.Ta(PCRproduct1, prtm[pn[i]], prtm[pn[j]])); //  adding by sorting and no repeatitions
                                    nampls++;
                                }
                                if (PCRproduct2 > mnPCRsize && PCRproduct2 < mxPCRsize) {
                                    map.put(PCRproduct2, Melting.Ta(PCRproduct2, prtm[pn[i]], prtm[pn[j]]));
                                    nampls++;
                                }
                            } else {
                                if (op[i] == 0) {
                                    if (PCRmatch_alignment) {
                                        sr.append("_____________________________\n");
                                        PrintPrimerAlignment(op[i], n, prl[pn[i]], x5[i], x3[i], pn[i]);
                                        PrintPrimerAlignment(op[j], n, prl[pn[j]], x5[j], x3[j], pn[j]);
                                    }
                                    sr.append("\n>").append(pname[pn[i]]).append("\t").append(s1).append("\n").append("5'-").append(primers[pn[i]]);
                                    sr.append("\n>").append(pname[pn[j]]).append("\t").append(s2).append("\n").append("5'-").append(primers[pn[j]]).append("\n");
                                } else {
                                    if (PCRmatch_alignment) {
                                        sr.append("_____________________________\n");
                                        PrintPrimerAlignment(op[j], n, prl[pn[j]], x5[j], x3[j], pn[j]);
                                        PrintPrimerAlignment(op[i], n, prl[pn[i]], x5[i], x3[i], pn[i]);
                                    }
                                    sr.append("\n>").append(pname[pn[j]]).append("\t").append(s2).append("\n").append("5'-").append(primers[pn[j]]);
                                    sr.append("\n>").append(pname[pn[i]]).append("\t").append(s1).append("\n").append("5'-").append(primers[pn[i]]).append("\n");
                                }
                                if (PCRproduct1 > mnPCRsize && PCRproduct1 < mxPCRsize) {
                                    sr.append("Amplicon size: ").append(PCRproduct1).append("bp\t").append("Ta=").append(tools.NumberToSeq(Melting.Ta(PCRproduct1, prtm[pn[i]], prtm[pn[j]]), 0)).append("°C\n");
                                }
                                if (PCRproduct2 > mnPCRsize && PCRproduct2 < mxPCRsize) {
                                    sr.append("Amplicon size: ").append(PCRproduct2).append("bp\t").append("Ta=").append(tools.NumberToSeq(Melting.Ta(PCRproduct2, prtm[pn[i]], prtm[pn[j]]), 0)).append("°C\n");
                                }
                                sr.append("\n");
                            }
                        }
                    }
                }
            }

            if (ShowOnlyAmplicons && !map.isEmpty()) {
                java.util.Set set = map.entrySet();
                java.util.Iterator iterator = set.iterator();
                sr.append("\nTotal number of amplicons = ").append(nampls).append("\n");
                sr.append("Number of amplicons of unique size = ").append(map.size()).append(" :\n");
                while (iterator.hasNext()) {
                    java.util.Map.Entry me2 = (java.util.Map.Entry) iterator.next();
                    sr.append(me2.getKey()).append("bp\t").append("Ta=").append(me2.getValue()).append("°C\n");
                }
            }
        }
    }

    private void PrintResult2(LinkList lnk, int n) {
        int l = 0;
        int i = 0;
        int j = 0;
        int t = 0;
        int h = 0;
        int m = 0;
        int y = 0;
        int x = 0;
        int e = 2;
        int b = 0;
        int d = 0;
        int z = 0;

        int l1 = seq[n].length();
        int n0 = 0;
        int n1 = 0;
        int n2 = 0;
        int n3 = 0;
        int nx = -1;

        String s1 = "";
        String s2 = "";
        String s3 = "";
        char[] value;
        char v;
        int k = lnk.Amount();
        int[] op = lnk.getOriantation();
        int[] x3 = lnk.getX3();
        int[] x5 = lnk.getX5();
        int[] pn = lnk.getPrimers();

        if (!lnk.isEmpty()) {
            sr.append("__________________________________________________\n In silico Primer(s) search for: ").append(CurrentFileName).append("//").append(sname[n]).append(":\n\n");
        }

        //prn = skolko vsego elements in pattern       
        //prx = ot kakogo FASTA primera, proishodi etot primer                   
        //prxn = n elements in current pattern       
        //ptx1 = pattern location between x1       
        //ptx2 = pattern location to x2    
        int[] d1;
        int[] d2;
        for (z = 0; z < prn.length; z++) {
            n1 = 0;
            n3 = 0;
            h = prn[z];

            d1 = new int[h + 1];
            d2 = new int[h + 1];

            for (m = 1; m <= h; m++) {
                t = 0;
                d1[m] = -1;
                d2[m] = k;
                for (i = n1; i < k; i++) {
                    if (prx[pn[i]] > z) {
                        d2[m] = i;
                        break;
                    }
                    if (prx[pn[i]] == z) {
                        if (t == 0 && (prxn[pn[i]] == m || h == 1)) { // index element in pattern, 1... prn[p]   
                            d1[m] = i;
                            t = 1;
                        }
                        if (prxn[pn[i]] > m) { // index element in pattern, 1... prn[p]   
                            d2[m] = i;
                            break;
                        }
                    }
                }
                if (d1[m] == -1) {
                    n3 = -1;
                    break;
                }
                n1 = d1[m];
            }

            if (n3 == 0) {
                for (b = d1[1]; b < d2[1]; b++) {
                    int[] rxe = new int[h];
                    n0 = 0;
                    s1 = "";
                    rxe[0] = b;

                    for (m = 2; m <= h; m++) {
                        if (op[b] == 0) {
                            for (d = d1[m]; d < d2[m]; d++) {
                                if (op[d] == 0) {
                                    if ((x3[rxe[n0]] + ptx2[pn[rxe[n0]]] > x5[d]) && (x3[rxe[n0]] + ptx1[pn[rxe[n0]]] <= x5[d])) {
                                        rxe[++n0] = d;
                                        break;
                                    }
                                    if ((x3[rxe[n0]] + ptx2[pn[rxe[n0]]] < x5[d])) {
                                        break;
                                    }
                                } else {
                                    //   break;
                                }
                            }
                        }

                        if (op[b] == 1) {
                            for (d = d2[m] - 1; d >= d1[m]; d--) {
                                if (op[d] == 1) {// reverse
                                    if ((x3[rxe[n0]] < x5[d] + ptx2[pn[rxe[n0]]]) && (x3[rxe[n0]] >= x5[d] + ptx1[pn[rxe[n0]]])) {
                                        rxe[++n0] = d;
                                        break;
                                    }
                                    if ((x3[rxe[n0]] > x5[d] + ptx2[pn[rxe[n0]]])) {
                                        break;
                                    }
                                } else {
                                    //  break;
                                }
                            }
                        }
                    }
                    if (n0 == h - 1 && nx < z) {
                        nx = z;
                        sr.append(">").append(pname[prx[pn[b]]]).append("\n").append(oprimers[prx[pn[b]]]).append("\n\n");
                    }
                    if (n0 == h - 1 && op[b] == 0) {
                        x = x5[rxe[0]] - e;
                        y = x3[rxe[n0]] + e + 1;
                        if (x < 0) {
                            s1 = (ct > 0) ? seq[n].substring(-x) : tools.Strings(-x, ' ');
                            x = 0;
                        }
                        if (y > l1) {
                            y = l1;
                        }
                        s1 = s1 + seq[n].substring(x, y);
                        l = s1.length();

                        s2 = pr[pn[rxe[0]]];
                        for (j = 1; j <= n0; j++) {
                            m = x5[rxe[j]] - x3[rxe[j - 1]] - 1;
                            s2 = (m > 0) ? s2 + tools.Strings(m, ' ') + pr[pn[rxe[j]]] : s2 + pr[pn[rxe[j]]];
                        }
                        s2 = "5-" + s2 + "->";

                        value = new char[l];
                        value[0] = ' ';
                        value[1] = ' ';
                        n1 = 0;
                        n2 = 0;
                        n3 = 0;
                        for (t = 2; t < l; t++) {
                            value[t] = ' ';
                            v = s2.charAt(t);
                            if (v != ' ') {
                                n3 = tables.aDNA[v][tables.cdna[s1.charAt(t)]];
                                if (n3 > tables.smmin) {
                                    value[t] = (n3 > tables.smth) ? '|' : ':';
                                }
                                n2++;
                                n1 = n1 + n3;
                            }
                        }
                        n1 = n1 / n2;
                        sr.append("Position: ").append(Xconvertor(x5[rxe[0]], x3[rxe[n0]], l1, ct, 0)).append("\t").append(n1).append("%\n");
                        sr.append(s2).append("\n").append(value).append("\n").append(s1).append("\n\n");
                        nsites++;
                    }

                    if (n0 == h - 1 && op[b] == 1) { // Reverse direction
                        x = x3[rxe[n0]] - e;
                        y = x5[rxe[0]] + e + 1;
                        s3 = "";
                        if (x < 0) {
                            s1 = (ct > 0) ? seq[n].substring(-x) : tools.Strings(-x, ' ');
                            x = 0;
                        }
                        if (y > l1) {
                            s3 = tools.Strings(y - l1, ' ');
                            y = l1;
                        }
                        s1 = s3 + dna.ComplementDNA(s1 + seq[n].substring(x, y));
                        l = s1.length();

                        s2 = pr[pn[rxe[0]]];
                        for (j = 1; j <= n0; j++) {
                            m = x3[rxe[j - 1]] - x5[rxe[j]] - 1;
                            s2 = (m > 0) ? s2 + tools.Strings(m, ' ') + pr[pn[rxe[j]]] : s2 + pr[pn[rxe[j]]];
                        }
                        s2 = "5-" + s2 + "->";

                        value = new char[l];
                        value[0] = ' ';
                        value[1] = ' ';
                        n1 = 0;
                        n2 = 0;
                        n3 = 0;
                        for (t = 2; t < l; t++) {
                            value[t] = ' ';
                            v = s2.charAt(t);
                            if (v != ' ') {
                                n3 = tables.aDNA[v][tables.cdna[s1.charAt(t)]];
                                if (n3 > tables.smmin) {
                                    value[t] = (n3 > tables.smth) ? '|' : ':';
                                }
                                n2++;
                                n1 = n1 + n3;
                            }
                        }
                        n1 = n1 / n2;
                        sr.append("Position: ").append(Xconvertor(x5[rxe[0]], x3[rxe[n0]], l1, ct, 2)).append("\t").append(n1).append("%\n");
                        sr.append(s2).append("\n").append(value).append("\n").append(s1).append("\n\n");
                        nsites++;
                    }
                }
            }
        }

    }

    private class LinkList {

        private int[] t;
        private int[] op;
        private int[] x3;
        private int[] x5;
        private int[] pn;
        private int n;
        private int p;
        private int e;

        public LinkList() {
            e = 1;// extending size
            p = 0;
            n = e;
            op = new int[e];
            x3 = new int[e];
            x5 = new int[e];
            pn = new int[e];
        }

        public boolean isEmpty() {
            return p <= 0;
        }

        public int Amount() {
            return p;
        }

        public int[] getOriantation() {
            return op;
        }

        public int[] getX3() {
            return x3;
        }

        public int[] getX5() {
            return x5;
        }

        public int[] getPrimers() {
            return pn;
        }

        public void insert(int _x3, int _x5, int _prn, int _op) {
            if (p == n) {
                n = n + e;
                t = new int[n];
                System.arraycopy(x3, 0, t, 0, p);// arraycopy(Object src, int srcPos, Object dest, int destPos, int length)   
                x3 = t;
                t = new int[n];
                System.arraycopy(x5, 0, t, 0, p);
                x5 = t;
                t = new int[n];
                System.arraycopy(pn, 0, t, 0, p);
                pn = t;
                t = new int[n];
                System.arraycopy(op, 0, t, 0, p);
                op = t;
            }
            op[p] = _op; // orientation
            x3[p] = _x3;
            x5[p] = _x5;
            pn[p] = _prn;// primer number
            p++;
        }

        public void shellSort() {
            if (p < 2 & nprimers > 0) {
                return;
            }
//http://thilinasameera.wordpress.com/2011/06/01/sorting-algorithms-sample-codes-on-java-c-and-matlab/        
//http://en.wikipedia.org/wiki/Sorting_algorithm            
            int n0 = 0;
            int n1 = 0;
            int n2 = 0;
            int n3 = 0;

            int j = 0;
            int t = 0;
            int l = p;
            int f = l / 2;
            while (f > 0) {
                for (int i = f; i < l; i++) {
                    n0 = x3[i];
                    n1 = x5[i];
                    n2 = pn[i];
                    n3 = op[i];
                    j = i;
                    while ((j >= f && pn[j - f] > n2) || ((j >= f && pn[j - f] == n2 && x5[j - f] > n1) || (j >= f && pn[j - f] == n2 && op[j - f] > n3))) { // sorting by primer n & x5
                        pn[j] = pn[j - f];
                        x3[j] = x3[j - f];
                        x5[j] = x5[j - f];
                        pn[j] = pn[j - f];
                        op[j] = op[j - f];
                        j = j - f;
                    }
                    x3[j] = n0;
                    x5[j] = n1;
                    pn[j] = n2;
                    op[j] = n3;
                }
                f = (f / 2);
            }
        }
    }

    public boolean isNumeric(String s) {
        return s.matches("[-+]?\\d*\\.?\\d+");
    }
    private HashMap<String, Integer> px;
    private int[] pvx1;
    private int[] pvx2;

    private String pr[]; //primers  
    private String apr[];//antisense primers    
    private String rpr[];//reverse primers    
    private String cpr[];//complement primers    
    private int prl[];   //primers length   
    private int pa[];    //n-bases adding to primer         
    private int prtm[];  //primer Tm    
    private int prn[];   //n elements in pattern       
    private int prx[];   //primer index for pattern    
    private int prxn[];  //n elements in current pattern       
    private int ptx1[];  //pattern x1       
    private int ptx2[];  //pattern x2      
    private int prnk[];  //F-R primers pairs control      
    private int u1[];    //coordinates x, -> px(,,,,1)    
    private int u2[];    //index sequence for u1 
    private int k1[];    //n elements per kmer
    private int p1[];    //primers  index
    private String CurrentFileName = "//";
    public int nprimers;
    public int nseq;
    private String[] seq;
    private String[] sname;
    private String primers[];
    private String oprimers[];
    private String pname[];
    private boolean pcr_predict;
    private boolean alignment;
    private boolean PCRmatch_alignment;
    private boolean CalculatePCRproduct;
    private boolean SequenceExtract;
    private boolean ShowOnlyAmplicons;
    private boolean patsearch;
    private boolean bis;
    private boolean probe;
    private boolean RFpairs;
    private int mxPCRsize;
    private int mnPCRsize;
    private int ct;
    private int er3;
    private int kmer;
    private int hashkmer;
    private int nsites;
    private int tprl; // total primer length
    private StringBuilder sr;
}
