// Algorithm 7 (Huff curves); needs core.m, montgomery.m
// Generated from the reference implementation; names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

function HuffStep(a, b, M1, M2, N)
    h := (N - 1) div 2;
    // Line 1
    An := a^2 + b^2;
    C := a*b;
    // Line 2
    t1 := (a - b)*(M1 + M2);
    t2 := (a + b)*(M1 - M2);
    // Line 3
    Anp, Cp, Xp, Zp, st := MontgomeryStep(An, C, t1 + t2, t1 - t2, N);
    X := st[1];
    Z := st[2];
    // Line 4
    wn := (a + b)^N*(&*[X[i] + Z[i] : i in [1..h]])^4;
    wd := (a - b)^N*(&*[X[i] - Z[i] : i in [1..h]])^4;
    // Line 5
    ap := wn + wd;
    bp := wn - wd;
    // Line 6
    t1 := (ap + bp)*(Xp + Zp);
    t2 := (ap - bp)*(Xp - Zp);
    return ap, bp, t1 + t2, t1 - t2, [* a, b, ap, bp, st *];
end function;

function HuffPush(M1, M2, st, N)
    a := st[1];
    b := st[2];
    ap := st[3];
    bp := st[4];
    stm := st[5];
    // Line 2
    t1 := (a - b)*(M1 + M2);
    t2 := (a + b)*(M1 - M2);
    // Line 3
    X, Z := MontgomeryPush(t1 + t2, t1 - t2, stm, N);
    // Line 6
    t1 := (ap + bp)*(X + Z);
    t2 := (ap - bp)*(X - Z);
    return t1 + t2, t1 - t2;
end function;

function HuffCodomainOnly(a, b, M1, M2, N)
    h := (N - 1) div 2;
    // Codomain-only variant of Algorithm 7: line 3 only computes the multiples X_i, Z_i.
    // Line 1
    An := a^2 + b^2;
    C := a*b;
    // Line 2
    t1 := (a - b)*(M1 + M2);
    t2 := (a + b)*(M1 - M2);
    // Line 3
    Fn, Gn, X, Z, Delta := MontMultiples(An, C, t1 + t2, t1 - t2, N, false);
    // Line 4
    wn := (a + b)^N*(&*[X[i] + Z[i] : i in [1..h]])^4;
    wd := (a - b)^N*(&*[X[i] - Z[i] : i in [1..h]])^4;
    // Line 5
    return wn + wd, wn - wd;
end function;
