import java.text.NumberFormat;
import java.util.HashMap;

public final class InSilicoPCR {

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

    public void SetPrimerStat(boolean a) {
        primer_statistic = a;
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

    // Extra bases reported on both sides of every extracted amplicon (0 = amplicon only).
    public void SetSequenceExtractFlanks(int v) {
        extractflanks = (v < 0) ? 0 : Math.min(v, 50000);
    }

    public void SetProductMaxLength(int v) {
        // +1 so the strict 'product < mxPCRsize' test includes an amplicon of exactly v,
        // mirroring the inclusive lower bound in SetProductMinLength.
        mxPCRsize = (v < 50000) ? v + 1 : 50001;
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
        kmer = 12;
        er3 = err3;

        if (er3 < 1) {
            er3 = 1;
        }
        if (er3 > 10) {
            er3 = 10;
        }

        if (er3 > 1) {
            kmer = 9;
        }

        ct = is_circle ? kmer : 0;
        if (!patsearch) {
            LoadPrimers1();
            if (probe) {
                PreparedPrimers(-1);
            } else {
                PreparedPrimers(16);
            }
        } else {
            LoadPrimers2(oriprimers);
            PreparedPrimers(-1);
        }
    }

    public void Run() {
        sr = new StringBuilder(1000);
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
        pstat = new int[nprimers];
        pstatseq = new int[nprimers];

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
        int k = 0;
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
            if (h > 0) {
                for (int i = 0; i < h; i++) {
                    if (!isNumeric(s[i])) {
                        if (tools.DNAtest(s[i]) == 1) {
                            s[i] = dna.DNA(s[i]);
                            l = s[i].length();
                            if (l > 11) {
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
                        s[0] = dna.DNA(s[0]);
                        l = s[0].length();
                        if (l > 11) {
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
            String[] s = prms[j].split("[\\[\\]\\-–]"); //String[] s = prms[j].split("\\p{Punct}");//http://docs.oracle.com/javase/6/docs/api/java/util/regex/Pattern.html#sum
            //String[] s = prms[j].split("(?<=\\p{Punct})");//http://docs.oracle.com/javase/6/docs/api/java/util/regex/Pattern.html#sum  "(?<=\\p{Punct})"
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
                            if (s[i].contains("@")) {
                                s[i] = dna.ComplementDNA(dna.DNA(s[i]));
                            } else {
                                s[i] = dna.DNA(s[i]);
                            }
                            l = s[i].length();
                            if (l > 11) {
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
        pstat = new int[nprimers];
        pstatseq = new int[nprimers];

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
        BuildCodeIndex();
    }

    // 2-bit codes for the four plain bases; -1 for everything else, including the degenerate
    // codes and the LNA/inosine spellings. Deliberately NOT folded: mapping u onto t or i onto g
    // here would make a primer and a target match where the strings do not.
    private static final byte[] B2 = new byte[128];

    static {
        for (int c = 0; c < 128; c++) {
            B2[c] = -1;
        }
        B2['a'] = 0;
        B2['c'] = 1;
        B2['g'] = 2;
        B2['t'] = 3;
    }

    // Builds the integer-keyed index over exactly those k-mers of px that are plain acgt. The scan
    // uses it only for windows it has proved are plain acgt too, so the two agree key for key; any
    // other window still goes through px by string, which costs an allocation but is rare.
    private void BuildCodeIndex() {
        pxc = new LongIntHashMap(Math.max(16, px.size() * 2));
        for (java.util.Map.Entry<String, Integer> e : px.entrySet()) {
            String k = e.getKey();
            if (k.length() != kmer) {
                continue;
            }
            int code = 0;
            boolean plain = true;
            for (int t = 0; t < kmer; t++) {
                char c = k.charAt(t);
                int v = (c < 128) ? B2[c] : -1;
                if (v < 0) {
                    plain = false;
                    break;
                }
                code = (code << 2) | v;
            }
            if (plain) {
                pxc.put(code, e.getValue());
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
            int v = l2 < 22 ? l2 - kmer - 1 : 20 - kmer; // test  
            v = v - Err3;

            if (y > -1) // direct Forward
            {
                x3 = x - y + l2 - 1;
                x5 = x3 - l2 + 1;
                r = 0;
                l = 0;
                if (x3 > ct - 1) {

                    // The old guard 'pvx1[d] < x3' assumed one primer's detections arrive in
                    // increasing x3. They do not: x3 is invariant across the k-mer offsets of a
                    // single site, but sites are visited in x = s + y order, so a genuine site
                    // found at a large offset after a nearer one was discarded outright, before
                    // any alignment ran. Repeats are collapsed after sorting instead, in dedup(),
                    // where equal (pn, x5, op) triples are adjacent.
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
                        if (tables.aDNA[c1][tables.cdna[c2]] < tables.sim24) { //<= tables.smth
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

                    if (probe || (!probe && r <= Err3)) {
                        e = 1;
                        // Mirrors the reverse branch: 'y - t > 0' stopped at t = y-1, so the
                        // primer's 5'-most base was never compared and a mismatch there went
                        // uncounted, while the same duplex written on the other strand did count
                        // it. 'x - t > 0' likewise hid base 0 of the target.
                        for (t = 1; t <= y; t++) {
                            if (x - t < 0) {
                                break;
                            }
                            c1 = pr[d].charAt(y - t);
                            c2 = seq[j].charAt(x - t);
                            if (tables.aDNA[c1][tables.cdna[c2]] < tables.sim24) { //<= tables.smmax
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
                        if (l > v) {
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

                    // Same as the forward branch: see the note on pvx1 above.
                    e = 1;
                    for (t = 1; t < y + 1; t++) {
                        if (x - t < 0) {
                            break;
                        }
                        c1 = rpr[d].charAt(y - t);
                        c2 = seq[j].charAt(x - t);
                        if (tables.aDNA[c1][c2] < tables.sim24) {
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
                    if (probe || (!probe && r <= Err3)) {
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
                            if (tables.aDNA[c1][c2] < tables.sim24) {
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
                        if (l > v) {
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
        int[] cb = new int[kmer];   // 2-bit value of each position in the current window, -1 if not acgt

        String[] R1;
        String s = "";
        for (int j = 0; j < nseq; j++) {
            nsites = 0;
            if (seq[j] != null) {
                l1 = seq[j].length();
                if (l1 > kmer) {

                    LinkList xr = new LinkList();

                    na = 1;
                    // cb/code/bad track the same window as ax/na, but in 2-bit form. They are
                    // maintained at the top of the loop, before any 'continue': the ax shift below
                    // is skipped on a miss (harmlessly, since na == 1 means every ax is 1), and a
                    // rolling code cannot be skipped that way without falling out of step.
                    int code = 0;
                    int bad = 0;
                    int slot = kmer - 1;
                    final int mask = (kmer >= 16) ? -1 : ((1 << (kmer + kmer)) - 1);
                    for (t = 0; t < kmer - 1; t++) {
                        char c0 = seq[j].charAt(t);
                        ax[t] = tables.t3[c0][0];
                        na = na * ax[t];
                        int v0 = (c0 < 128) ? B2[c0] : -1;
                        cb[t] = v0;
                        if (v0 < 0) {
                            bad++;
                            v0 = 0;
                        }
                        code = ((code << 2) | v0) & mask;
                    }
                    // A circular template needs exactly kmer-1 extra windows (the ones straddling
                    // the origin). The old bound l1 + ct ran one further, to x == l1, where the
                    // wrap branch below built a k-mer from a different offset than ax[]/na describe
                    // and AddRecord was called past the end of the molecule. Written this way, not
                    // as l1 + ct - 1, because ct == 0 must keep the bound at l1: shortening it there
                    // would silently drop the last window of every linear sequence.
                    for (int i = kmer - 1; i < l1 + (ct > 0 ? ct - 1 : 0); i++) {
                        x = 1 + i - kmer;
                        boolean wrapped = i > l1 - 1; // circle, ct=kmer
                        char cin = wrapped ? seq[j].charAt(i - l1) : seq[j].charAt(i);
                        ax[kmer - 1] = tables.t3[cin][0];
                        na = na * ax[kmer - 1];

                        if (i >= kmer && cb[slot] < 0) { // the base dropping out of the window
                            bad--;
                        }
                        int vin = (cin < 128) ? B2[cin] : -1;
                        cb[slot] = vin;
                        if (vin < 0) {
                            bad++;
                            vin = 0;
                        }
                        code = ((code << 2) | vin) & mask;
                        if (++slot == kmer) {
                            slot = 0;
                        }
                        if (na < dna.maxdna) {
                            if (na == 1 && bad == 0) {
                                // The whole window is plain acgt, so the rolling code identifies it
                                // exactly and no k-mer string is built at all - this is the path
                                // every base of an ordinary target takes.
                                int p0 = pxc.get(code);
                                if (p0 < 0) { // w starts at -1 and is incremented before put, so 0 is a
                                    continue; // real index; only NO_VALUE (Integer.MIN_VALUE) means absent
                                }
                                AddRecord(x, j, p0, l1, Err3, xr);
                                na = na / ax[0];
                                for (t = 0; t < kmer - 1; t++) {
                                    ax[t] = ax[t + 1];
                                }
                                continue;
                            }
                            // Anything else still needs the k-mer as text: a window that is not
                            // degenerate but is not plain acgt either, or one to be expanded.
                            if (wrapped) {
                                // x < l1 always holds now, so only the joined form remains
                                s = seq[j].substring(x) + seq[j].substring(0, i - l1 + 1);
                            } else {
                                s = seq[j].substring(x, x + kmer);
                            }
                            if (na == 1) {
                                Integer p0 = px.get(s); // one probe: containsKey + get walked the bucket twice
                                if (p0 == null) {
                                    continue;
                                }
                                p = p0;
                                AddRecord(x, j, p, l1, Err3, xr);
                            } else {
                                if ((ax[kmer - 1] * ax[kmer - 2] * ax[kmer - 3]) < 64 && (ax[0] * ax[1] * ax[2]) < 64) {
                                    R1 = dna.ConvertDNA(s); // use the (possibly circular-wrapped) k-mer 's'; raw substring would read past the junction
                                    for (String R11 : R1) {
                                        Integer p0 = px.get(R11);
                                        if (p0 == null) {
                                            continue;
                                        }
                                        p = p0;
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
                            xr.dedup();
                            PrintResult2(xr, j);
                            sr.append("\nAmount of sites: ").append(nsites).append("\n\n");

                        } else {
                            xr.shellSort();
                            xr.dedup();
                            PrintResult1(xr, j);
                        }
                    }
                }
            }
        }

        //primers statistic for the current file (summed over all its sequences)
        if (primer_statistic) {
            sr.append("\nPrimerID\tSequence\tHits\n");
            Integer[] order = new Integer[nprimers > 1 ? nprimers - 1 : 0];
            for (int j = 1; j < nprimers; j++) {
                order[j - 1] = j;
            }
            java.util.Arrays.sort(order, (a, b) -> Integer.compare(pstat[b], pstat[a])); // by Hits, descending (stable for ties)
            if (!patsearch) {
                for (int j : order) {
                    sr.append(pname[j]).append("\t").append(pr[j]).append("\t").append(pstat[j]).append("\n");
                }
            } else {
                for (int j : order) {
                    sr.append(pname[prx[j]]).append("\t").append(pr[j]).append("\t").append(pstat[j]).append("\n");
                }
            }
            sr.append("\n");
            // The global summary is fed by the caller via publishGlobalStat(), not from here:
            // with several files in flight at once this method runs on a worker thread.
        }
    }

    // ---- Global (cross-file) primer occurrence statistics ----
    // When many files are analysed in one run, each file already prints its own
    // per-primer table. These members additionally sum the frequencies per primer
    // across every file so a single consolidated table can be printed at the very end.
    private static java.util.LinkedHashMap<String, GStat> globalStat;
    private static long globalSeqTotal;
    private static long globalFileTotal;

    private static final class GStat {

        String id;   // primer ID
        String seq;  // primer sequence (as shown in the per-file table)
        long hits;   // total binding sites over all files
        long seqs;   // number of sequences (over all files) the primer was found in
        long files;  // number of files the primer was found in
    }

    // Start collecting a cross-file summary. Call once before processing a set of files.
    public static void beginGlobalStat() {
        globalStat = new java.util.LinkedHashMap<>();
        globalSeqTotal = 0;
        globalFileTotal = 0;
    }

    // Render the consolidated, human-readable table (sorted by Hits, descending).
    // Returns "" if inactive or empty. Does not clear the collected data (beginGlobalStat resets it),
    // so it may be called together with getGlobalStatTsv().
    public static String getGlobalStatReport() {
        if (globalStat == null || globalStat.isEmpty()) {
            return "";
        }
        StringBuilder b = new StringBuilder(256);
        b.append("__________________________________________________\n");
        b.append("Global primer occurrence statistics (").append(globalFileTotal)
                .append(" file(s), ").append(globalSeqTotal).append(" sequence(s))\n");
        appendStatTable(b);
        b.append("\n");
        return b.toString();
    }

    // Render the consolidated table as plain TSV (header + rows, sorted by Hits, descending),
    // for a stand-alone *.stats.tsv file. Returns "" if inactive or empty.
    public static String getGlobalStatTsv() {
        if (globalStat == null || globalStat.isEmpty()) {
            return "";
        }
        StringBuilder b = new StringBuilder(256);
        appendStatTable(b);
        return b.toString();
    }

    // Append the shared header + one row per primer (sorted by Hits, descending) to b.
    private static void appendStatTable(StringBuilder b) {
        b.append("PrimerID\tSequence\tHits\tSequences\tFiles\n");
        for (GStat g : sortedGlobal()) {
            b.append(g.id).append("\t").append(g.seq).append("\t")
                    .append(g.hits).append("\t").append(g.seqs).append("\t").append(g.files).append("\n");
        }
    }

    private static java.util.List<GStat> sortedGlobal() {
        java.util.List<GStat> list = new java.util.ArrayList<>(globalStat.values());
        list.sort((a, b) -> Long.compare(b.hits, a.hits)); // by Hits, descending (stable for ties)
        return list;
    }

    // Add this instance's per-file statistics into the cross-file summary (no-op if inactive).
    // Folds this run's per-primer tallies into the cross-file summary. Call it from the thread
    // that owns the summary, once per processed file and in file order; never from a worker.
    public void publishGlobalStat() {
        if (globalStat == null || !primer_statistic) {
            return;
        }
        globalSeqTotal += nseq;
        globalFileTotal++;
        if (!patsearch) {
            for (int j = 1; j < nprimers; j++) {
                addGlobal(pname[j], pr[j], pstat[j], pstatseq[j]);
            }
        } else {
            for (int j = 1; j < nprimers; j++) {
                addGlobal(pname[prx[j]], pr[j], pstat[j], pstatseq[j]);
            }
        }
    }

    private static void addGlobal(String id, String seq, int hits, int seqs) {
        id = (id == null) ? "" : id.trim();
        String key = id + "\t" + seq;
        GStat g = globalStat.get(key);
        if (g == null) {
            g = new GStat();
            g.id = id;
            g.seq = seq;
            globalStat.put(key, g);
        }
        g.hits += hits;
        g.seqs += seqs;
        if (hits > 0) {
            g.files++;
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
                s1 = (ct > 0) ? seq[n].substring(l1 + x) : tools.Strings(-x, ' ');
                x = 0;
            }

            if (y > l1 - 1) {
                s2 = (ct > 0) ? seq[n].substring(0, 1 + y - l1) : tools.Strings(1 + y - l1, ' ');
                y = l1 - 1;
            }
            s = dna.AntisenseDNA(s1 + seq[n].substring(x, y + 1) + s2);
            tm = Melting.Tmelting(pr[pn], s, 0.2d, 55d, 1, 0, 4);

            value = new char[l + 2];
            value[0] = ' ';
            value[1] = ' ';

            n1 = 0;
            for (t = 0; t < l; t++) {
                value[t + 2] = ' ';
                n3 = tables.aDNA[pr[pn].charAt(t)][s.charAt(t)];
                if (n3 > tables.sim24) {
                    value[t + 2] = n3 > tables.smmin ? '|' : ':';
                }
                n1 = n1 + n3;
            }
            n1 = n1 / l;

            x = x5 - 2;
            y = x3 + 4;
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
            s = dna.AntisenseDNA(s1 + seq[n].substring(x, y + 1) + s2);
            s1 = (pa[pn] > 0) ? tools.Strings(pa[pn], ' ') + "5-" : "5-";
            sr.append("Position: ").append(Xconvertor(x5, x3, l1, ct, 0)).append("\t").append(n1).append("%\tTm=").append(tools.NumberToSeq(tm, 1)).append("\n");
            sr.append(s1).append(primers[pn]).append("->\n").append(new String(value)).append("\n").append(s).append("\n");
        } else 
            ///////////////// reverse primer
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
                if (n3 > tables.sim24) {
                    value[t + 2] = n3 > tables.smmin ? '|' : ':';
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
            sr.append("Position: ").append(Xconvertor(x5, x3, l1, ct, 1)).append("\t").append(n1).append("%\tTm=").append(tools.NumberToSeq(tm, 1)).append("\n");
            sr.append("<-").append(rpr[pn]).append("-5\n").append(new String(value)).append("\n").append(s).append("\n");
        }
    }

    private void PrintResult1(LinkList lnk, int n) {
        NumberFormat nf = NumberFormat.getInstance();
        int i = 0;
        int j = 0;
        int x1 = 0;
        int x2 = 0;
        int nampls = 0;

        int k = lnk.Amount();
        int l1 = seq[n].length();
        int n0 = -1;

        String s1 = "";
        String s2 = "";

        int[] op = lnk.getOriantation();
        int[] x3 = lnk.getX3();
        int[] x5 = lnk.getX5();
        int[] pn = lnk.getPrimers();

        if (primer_statistic) {
            int last = -1;
            for (j = 0; j < k; j++) {
                pstat[pn[j]]++;
                if (pn[j] != last) { // hits are grouped by primer (sorted): count each primer once per sequence
                    pstatseq[pn[j]]++;
                    last = pn[j];
                }
            }
        }

        java.util.Map<Integer, Integer> map = new java.util.TreeMap<>(java.util.Collections.reverseOrder());
        java.util.Map<Integer, String> mapNames = new java.util.TreeMap<>(java.util.Collections.reverseOrder());
        java.util.Set<Long> seen2 = new java.util.HashSet<>(); // amplicon regions already extracted

        sr.append("__________________________________________________\nIn silico PCR, primers search for: ").append(CurrentFileName).append("\n").append(sname[n]).append("\n");
        sr.append("Target sequence length = ").append(nf.format(l1)).append(" nt").append("\n\n");

        if (alignment) {
            for (i = 0; i < k; i++) {
                if (n0 != pn[i]) {
                    n0 = pn[i];
                    sr.append(pname[pn[i]]).append("\t").append(primers[pn[i]]).append("\n\n");
                }
                PrintPrimerAlignment(op[i], n, prl[pn[i]], x5[i], x3[i], pn[i]);
                sr.append("\n");
            }
        } else {
            for (j = 0; j < k; j++) {
                if (n0 != pn[j]) {
                    sr.append("\n").append(pname[pn[j]]);
                    n0 = pn[j];
                }
                if (op[j] == 0) {
                    sr.append("\t").append(nf.format(x5[j])).append("->").append(nf.format(x3[j]));
                } else {
                    sr.append("\t").append(nf.format(x3[j])).append("<-").append(nf.format(x5[j]));
                }
            }
            sr.append("\n");
        }

        if (pcr_predict || SequenceExtract || ShowOnlyAmplicons) {
            // A pair can only be reported when its product length falls inside
            // (mnPCRsize, mxPCRsize), and on a linear template that length is a plain difference
            // of the two 5' positions. So the partners of site i are exactly the sites whose x5
            // lies in one interval, and they can be looked up instead of found by scanning every
            // one of the k*(k-1)/2 combinations. A circular template also pairs across the origin,
            // which is not one interval, so it keeps the exhaustive loop - as does a small k,
            // where sorting would cost more than the scan it saves.
            // Order is preserved exactly: candidates are filtered to j > i and walked ascending,
            // so pairs are still visited in (i, j) order and the report is unchanged.
            boolean windowed = (ct == 0) && (k > 128);
            long[] ord = null;   // ((long) x5 << 32) | siteIndex, sorted -> ordered by x5
            int[] cj = null;     // candidate partners of the current i, ascending
            if (windowed) {
                ord = new long[k];
                for (int t = 0; t < k; t++) {
                    ord[t] = ((long) x5[t] << 32) | t;
                }
                java.util.Arrays.sort(ord);
                cj = new int[64];
            }
            for (i = 0; i < k - 1; i++) {
                int cn;
                if (windowed) {
                    long lo;    // x5[j] must be strictly inside (lo, hi)
                    long hi;
                    if (op[i] == 0) {   // product = 1 + x5[j] - x5[i]
                        lo = (long) x5[i] + mnPCRsize - 1;
                        hi = (long) x5[i] + mxPCRsize - 1;
                    } else {            // product = 1 + x5[i] - x5[j]
                        lo = (long) x5[i] - mxPCRsize + 1;
                        hi = (long) x5[i] - mnPCRsize + 1;
                    }
                    cn = 0;
                    for (int q = lowerBound(ord, lo + 1); q < ord.length && (ord[q] >> 32) <= hi - 1; q++) {
                        int t = (int) ord[q];
                        if (t > i) {
                            if (cn == cj.length) {
                                cj = java.util.Arrays.copyOf(cj, cn << 1);
                            }
                            cj[cn++] = t;
                        }
                    }
                    java.util.Arrays.sort(cj, 0, cn);   // back into site order, as the plain loop had them
                } else {
                    cn = k - i - 1;
                }
                for (int c = 0; c < cn; c++) {
                    j = windowed ? cj[c] : (i + 1 + c);
                    if ((op[i] != op[j]) && prnk[pn[i]] == prnk[pn[j]]) {
                        int PCRproduct1 = 0;
                        int PCRproduct2 = 0;
                        int x0 = 0;
                        if (op[i] == 0) { //  Fi-> <-jR
                            x0 = 1;
                            x1 = x5[i];
                            x2 = x5[j];
                            PCRproduct1 = 1 + x5[j] - x5[i]; //<-jR Fi->  F(i)->L1___1<-(j) ' Fi-> <-jR   
                            if (ct > 0 && x5[i] > x5[j]) {
                                // Inclusive, exactly like PCRproduct1: the reverse site's x5 is its
                                // rightmost base, so the wrapping product spans 1 + (l1 - x5[i]) + x5[j].
                                // No else-branch: for one F/R pair on a circle there is a single
                                // product, and the old else-if only re-stated this expression for the
                                // case where the direct product already existed.
                                PCRproduct2 = 1 + l1 - x5[i] + x5[j];
                            }
                        } else {  // Fj-> <-iR
                            x0 = 2;
                            x1 = x5[j];
                            x2 = x5[i];
                            PCRproduct1 = 1 + x5[i] - x5[j]; //<-iR Fj->  F(j)->L1---1<-(i) ' Fj-> <-iR                    
                            if (ct > 0 && x5[j] > x5[i]) {
                                PCRproduct2 = 1 + l1 - x5[j] + x5[i];
                            }
                        }

                        if (ct > 0 && PCRproduct1 > l1) {
                            // Both sites straddle the origin: the forward one was recorded with a
                            // negative x5 and the reverse one with an x5 past l1, so the plain
                            // difference above counted one extra turn of the molecule. A product
                            // cannot be longer than the circle it sits on.
                            PCRproduct1 = PCRproduct1 - l1;
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
                            // On a circle x5 is legitimately negative (a site straddling the origin)
                            // or past l1 (a reverse site straddling it, since x5 = x - y + l2 - 1);
                            // clamping those to 0 or l1 silently shortened the extracted record while
                            // the header kept the full size. Only the start is normalised here - the
                            // end and the wrap flag are derived from the start and the product length
                            // in ExtractAmplicon, because the two coordinates alone no longer say
                            // which way round the molecule the amplicon runs.
                            if (ct > 0) {
                                x1 = Mod(x1, l1);
                            } else {
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
                            }

                            if (SequenceExtract) {
                                if (PCRproduct1 > mnPCRsize && PCRproduct1 < mxPCRsize) {
                                    if (ct > 0) {
                                        ExtractAmplicon(n, x1, PCRproduct1, pn[i], pn[j], nf, seen2);
                                    } else if (seen2.add(RegionKey(x1, x2, false))) {
                                        AppendAmplicon(n, x1, x2, PCRproduct1, false, pn[i], pn[j], nf);
                                    }
                                }
                                // The origin-spanning product exists only when x5[i] > x5[j], and the
                                // direct product is then non-positive, so the two are mutually
                                // exclusive and cannot collide on the same region key.
                                if (PCRproduct2 > mnPCRsize && PCRproduct2 < mxPCRsize) {
                                    ExtractAmplicon(n, x1, PCRproduct2, pn[i], pn[j], nf, seen2);
                                }
                            }

                            if (ShowOnlyAmplicons) {
                                if (PCRproduct1 > mnPCRsize && PCRproduct1 < mxPCRsize) {
                                    map.put(PCRproduct1, Melting.Ta(PCRproduct1, prtm[pn[i]], prtm[pn[j]])); //  adding by sorting and no repeatitions
                                    mapNames.put(PCRproduct1, pname[pn[i]].trim() + " / " + pname[pn[j]].trim());
                                    nampls++;
                                }
                                if (PCRproduct2 > mnPCRsize && PCRproduct2 < mxPCRsize) {
                                    map.put(PCRproduct2, Melting.Ta(PCRproduct2, prtm[pn[i]], prtm[pn[j]]));
                                    mapNames.put(PCRproduct2, pname[pn[i]].trim() + " / " + pname[pn[j]].trim());
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
                                    sr.append("Amplicon size: ").append(PCRproduct1).append("bp\n");
                                }
                                if (PCRproduct2 > mnPCRsize && PCRproduct2 < mxPCRsize) {
                                    sr.append("Amplicon size: ").append(PCRproduct2).append("bp\n");
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
                sr.append("Number of amplicons of unique size = ").append(map.size()).append("\n");
                while (iterator.hasNext()) {
                    java.util.Map.Entry me2 = (java.util.Map.Entry) iterator.next();
                    sr.append(me2.getKey()).append("bp\t").append("Ta=").append(me2.getValue()).append("\t").append(mapNames.get(me2.getKey())).append("\n");
                }
            }
        }
    }

    // First position in the x5-ordered index whose x5 is not below x5min.
    private static int lowerBound(long[] ord, long x5min) {
        int lo = 0;
        int hi = ord.length;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if ((ord[mid] >> 32) < x5min) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        return lo;
    }

    // Writes one amplicon of a circular template, given where it starts and how long it is.
    // Those two are what the caller actually knows: the end coordinate and the wrap flag are
    // derived here so they cannot disagree with the size printed in the header. Deriving the flag
    // instead by comparing the two wrapped coordinates fails exactly when the product reaches l1 -
    // the residues then compare equal or the wrong way round, and the record comes out short.
    private void ExtractAmplicon(int n, int start, int size, int p1, int p2, NumberFormat nf, java.util.Set<Long> seen) {
        int l1 = seq[n].length();
        int take = Math.min(size, l1);   // a product cannot run more than once round the molecule
        int end = Mod(start + take - 1, l1);
        boolean wrap = start + take > l1;
        if (seen.add(RegionKey(start, end, wrap))) {
            AppendAmplicon(n, start, end, size, wrap, p1, p2, nf);
        }
    }

    // Key identifying one extracted region, so the same amplicon is written only once per sequence.
    private static long RegionKey(int x1, int x2, boolean wrap) {
        return (wrap ? 1L << 62 : 0L) | ((long) x1 << 31) | x2;
    }

    // Non-negative remainder: keeps a coordinate inside a circular template.
    private static int Mod(int v, int m) {
        int r = v % m;
        return (r < 0) ? r + m : r;
    }

    // `len` bases starting at the 0-based position `from`: a circular template wraps around
    // the origin, a linear one is simply clipped at its ends (a shorter flank is then returned).
    private static String SubSequence(String s, int from, int len, boolean circular) {
        int l1 = s.length();
        if (len <= 0 || l1 == 0) {
            return "";
        }
        if (!circular) {
            int a = Math.max(from, 0);
            int b = Math.min(from + len, l1);
            return (a < b) ? s.substring(a, b) : "";
        }
        if (len > l1) {
            len = l1;
        }
        int a = Mod(from, l1);
        int b = a + len;
        return (b <= l1) ? s.substring(a, b) : s.substring(a) + s.substring(0, b - l1);
    }

    // Writes one extracted amplicon in FASTA form. With extractflanks > 0 the record is extended
    // by up to that many bases on both sides (lower case), the amplicon itself staying upper case;
    // wrap = true marks a product of a circular template that runs across the origin.
    private void AppendAmplicon(int n, int x1, int x2, int size, boolean wrap, int p1, int p2, NumberFormat nf) {
        String s = seq[n];
        int l1 = s.length();
        boolean circular = ct > 0;

        // x2 is the amplicon's last base, so both forms are inclusive of it.
        String core = wrap ? s.substring(x1) + s.substring(0, Math.min(1 + x2, l1))
                : s.substring(x1, Math.min(1 + x2, l1));
        int endx = Math.min(1 + x2, l1); // first base after the amplicon, 0-based

        int f = extractflanks;
        if (circular) { // on a circle the flanks must not run into the amplicon from the other side
            f = Math.min(f, Math.max((l1 - core.length()) / 2, 0));
        }
        String left = SubSequence(s, x1 - f, f, circular);
        String right = SubSequence(s, endx, f, circular);

        int a1 = 1 + x1;                                // amplicon start/end, 1-based
        int a2 = Math.min(1 + x2, l1);                  // now that x2 is included, no special case

        // Kept from the previous version: the annealing temperature can be added to the header with
        // .append("\tTa=").append(tools.NumberToSeq(Melting.Ta(size, prtm[p1], prtm[p2]), 0)).append("°C")
        sr.append("\n");
        sr.append(pname[p1]).append("\t").append(primers[p1]).append("\n");
        sr.append(pname[p2]).append("\t").append(primers[p2]).append("\n");

        if (left.isEmpty() && right.isEmpty()) {
            sr.append(">").append(nf.format(a1)).append("-").append(nf.format(a2)).append("(").append(size).append(" bp)").append("\n");
            sr.append(core).append("\n\n");
        } else {
            int r1 = 1 + Mod(x1 - left.length(), l1);   // extracted region start/end, 1-based
            int r2 = 1 + Mod(endx + right.length() - 1, l1);
            sr.append(">").append(nf.format(r1)).append("-").append(nf.format(r2))
                    .append("(").append(left.length() + core.length() + right.length()).append(" bp)")
                    .append("\tamplicon: ").append(nf.format(a1)).append("-").append(nf.format(a2))
                    .append("(").append(size).append(" bp)")
                    .append("\tflanks: ").append(left.length()).append("+").append(right.length()).append(" bp\n");
            sr.append(left).append(core.toUpperCase(java.util.Locale.ROOT)).append(right).append("\n\n");
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

        if (primer_statistic) {
            int last = -1;
            for (j = 0; j < k; j++) {
                pstat[pn[j]]++;
                if (pn[j] != last) { // hits are grouped by primer (sorted): count each primer once per sequence
                    pstatseq[pn[j]]++;
                    last = pn[j];
                }
            }
        }

        if (!lnk.isEmpty()) {
            sr.append("__________________________________________________\nLinked search for: ").append(CurrentFileName).append("\n").append(sname[n]).append("\n");
            sr.append("Target sequence length = ").append(l1).append(" nt").append("\n\n");
        }

        //prn = skolko vsego elements in pattern       
        //prx = ot kakogo FASTA primera, proishodi etot primer                   
        //prxn = n elements in current pattern       
        //ptx1 = pattern location between x1       
        //ptx2 = pattern location to x2    
        int[] d1;
        int[] d2;
        for (z = 1; z < prn.length; z++) {
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
                        sr.append(pname[prx[pn[b]]]).append("\n").append(oprimers[prx[pn[b]]]).append("\n\n");
                    }
                    if (n0 == h - 1 && op[b] == 0) {
                        x = x5[rxe[0]] - e;
                        y = x3[rxe[n0]] + e + 1;
                        s1 = ""; // else the previous window's flank is prepended below
                        if (x < 0) {
                            s1 = (ct > 0) ? seq[n].substring(l1 + x) : tools.Strings(-x, ' ');
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
                        for (t = 2; t < l; t++) {
                            value[t] = ' ';
                            v = (t < s2.length()) ? s2.charAt(t) : ' '; // s2 may be shorter than the s1 window (circular/edge); avoid OOB
                            if (v != ' ') {
                                n3 = tables.aDNA[v][tables.cdna[s1.charAt(t)]];
                                if (n3 > tables.sim24) {
                                    value[t] = (n3 > tables.smmin) ? '|' : ':';
                                }
                                n2++;
                                n1 = n1 + n3;
                            }
                        }
                        n1 = (n2 > 0) ? n1 / n2 : 0; // guard against an all-gap window
                        sr.append("Position: ").append(Xconvertor(x5[rxe[0]], x3[rxe[n0]], l1, ct, 0)).append("\t").append(n1).append("%\n");
                        sr.append(s2).append("\n").append(value).append("\n").append(s1).append("\n\n");
                        nsites++;
                    }

                    if (n0 == h - 1 && op[b] == 1) { // Reverse direction
                        x = x3[rxe[n0]] - e;
                        y = x5[rxe[0]] + e + 1;
                        s1 = ""; // else the previous window's flank is prepended below
                        s3 = "";
                        if (x < 0) {
                            s1 = (ct > 0) ? seq[n].substring(l1 + x) : tools.Strings(-x, ' ');
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
                        for (t = 2; t < l; t++) {
                            value[t] = ' ';
                            v = (t < s2.length()) ? s2.charAt(t) : ' '; // s2 may be shorter than the s1 window (circular/edge); avoid OOB
                            if (v != ' ') {
                                n3 = tables.aDNA[v][tables.cdna[s1.charAt(t)]];
                                if (n3 > tables.sim24) {
                                    value[t] = (n3 > tables.smmin) ? '|' : ':';
                                }
                                n2++;
                                n1 = n1 + n3;
                            }
                        }
                        n1 = (n2 > 0) ? n1 / n2 : 0; // guard against an all-gap window
                        sr.append("Position: ").append(Xconvertor(x5[rxe[0]], x3[rxe[n0]], l1, ct, 2)).append("\t").append(n1).append("%\n");
                        sr.append(s2).append("\n").append(value).append("\n").append(s1).append("\n\n");
                        nsites++;
                    }
                }
            }
        }

    }

    private class LinkList {

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
                // Geometric growth: growing by e == 1 reallocated and copied all four arrays on
                // every accepted site, i.e. O(m^2) copying over m sites. The initial size stays 1,
                // a LinkList is built per sequence and an assembly may hold thousands of contigs.
                n = n + (n < 16 ? 16 : n >> 1);
                int[] t = new int[n];
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

        // Collapses repeated detections of the same site. Must run on the sorted arrays, where
        // equal (pn, x5, op) triples are adjacent; x3 is a function of those three, so it follows.
        public void dedup() {
            if (p < 2) {
                return;
            }
            int w = 1;
            for (int r = 1; r < p; r++) {
                if (pn[r] == pn[w - 1] && x5[r] == x5[w - 1] && op[r] == op[w - 1]) {
                    continue;
                }
                pn[w] = pn[r];
                x5[w] = x5[r];
                x3[w] = x3[r];
                op[w] = op[r];
                w++;
            }
            p = w;
        }

        public void shellSort() {
            if (p < 2 & nprimers > 0) {
                return;
            }
//http://thilinasameera.wordpress.com/2011/06/01/sorting-algorithms-sample-codes-on-java-c-and-matlab/        
//http://en.wikipedia.org/wiki/Sorting_algorithm            

            int l = p;
            int f = l / 2;
            while (f > 0) {
                for (int i = f; i < l; i++) {
                    int n0 = x3[i];
                    int n1 = x5[i];
                    int n2 = pn[i];
                    int n3 = op[i];
                    int j = i;
                    // Lexicographic by (pn, x5, op). The op clause used to be a separate
                    // disjunct, so it fired regardless of x5 and the three keys did not define a
                    // total order - the sort could leave equal-key neighbours apart, which any
                    // neighbour-based deduplication downstream depends on.
                    while (j >= f && (pn[j - f] > n2
                            || (pn[j - f] == n2 && (x5[j - f] > n1
                            || (x5[j - f] == n1 && op[j - f] > n3))))) { // sorting by primer n, then x5, then orientation
                        pn[j] = pn[j - f];
                        x3[j] = x3[j - f];
                        x5[j] = x5[j - f];
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
    private LongIntHashMap pxc; // the same k-mers keyed by their 2-bit code, for plain-acgt windows

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
    private int pstat[];
    private int pstatseq[]; //number of sequences (within the current file) in which the primer was found
    private boolean pcr_predict;
    private boolean alignment;
    private boolean PCRmatch_alignment;
    private boolean SequenceExtract;
    private int extractflanks; // flanking bases added to each extracted amplicon
    private boolean ShowOnlyAmplicons;
    private boolean patsearch;
    private boolean primer_statistic;
    private boolean probe;
    private boolean RFpairs;
    private int mxPCRsize;
    private int mnPCRsize;
    private int ct;
    private int er3;
    private int kmer;
    private int nsites;
    private int tprl; // total primer length
    private StringBuilder sr;
}
