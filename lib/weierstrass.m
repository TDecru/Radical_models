// Algorithm 3 (short Weierstrass) and its codomain-only variant; needs core.m
// Generated from the reference implementation; names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

function SWMultiples(a4, a6, X1, Z1, N, rs)
    // Lines 1-3: F_n, G_n and the shifted multiples (Delta_k, D_k). If rs = 1, line 2 rescales
    // (X_1 : Z_1) by F_n, so that D_2 = 1 (and F_n, G_n change accordingly).
    h := (N - 1) div 2;
    Delta := AssociativeArray();
    D := AssociativeArray();
    // Line 1
    t0 := Z1^2;
    t1 := X1^2;
    t2 := a4*t0;
    t3 := a6*t0*Z1;
    Fn := 4*(X1*(t1 + t2) + t3);
    Gn := 3*t1 + t2;
    // Line 2
    Delta[1] := 0;
    D[1] := 1;
    Delta[2] := (t1 - t2)^2 - X1*(Fn + 8*t3);
    D[2] := Fn;
    Xo := X1;
    Zo := Z1;
    if rs eq 1 then
        t5 := Fn^2;
        Xo := X1*Fn;
        Zo := Z1*Fn;
        Gn := Gn*t5;
        Fn := t5^2;
        D[2] := 1;
    end if;
    // Line 3
    for k in [2..h - 1] do
        t4 := Delta[k]^2;
        Delta[k + 1] := D[k - 1]*D[k]*(Fn*D[k] + 2*Gn*Delta[k]) - t4*Delta[k - 1];
        D[k + 1] := t4*D[k - 1];
    end for;
    return Fn, Gn, Delta, D, Xo, Zo;
end function;

function SWStep(a4, a6, X1, Z1, N)
    h := (N - 1) div 2;
    // Line 1
    Fn, Gn, Delta, D, X1s, Z1s := SWMultiples(a4, a6, X1, Z1, N, 1);
    // Line 4
    PD2 := Delta[h]^2*D[h - 1]*D[h];
    g, g2, g3, R, U, S := CoreRadical(Fn, 1, PD2, Delta, D, N);
    Rd, FnR, abar, bbar, cbar := CoreCubic(Fn, Gn, 1, Delta[2], D[2], g, g2, g3, R, U, S);
    // Line 5
    XV := X1s*Rd + FnR;
    ZV := Z1s*Rd;
    t1 := XV*ZV;
    t0 := ZV*abar;
    // Line 6
    t2 := t1*abar;
    t3 := t2^2;
    t4 := abar^2;
    a4p := t3*(9*bbar - 3*t4);
    a6p := t3*t2*(abar*(2*t4 - 9*bbar) + 27*cbar);
    xPp := t2*abar;
    st := [* X1s, Z1s, Delta, D, a4, a6, t0 *];
    return a4p, a6p, xPp, st;
end function;

function SWCommonDenominator(X1, Z1, Delta, D, N)
    // Kernel x-coordinates x_i = Xb[i]/Dc over a common denominator (running products of the D_k from both ends).
    // Line 7
    h := (N - 1) div 2;
    pre := AssociativeArray();
    suf := AssociativeArray();
    Xb := AssociativeArray();
    pre[1] := 1;
    for k in [2..h - 1] do
        pre[k] := pre[k - 1]*D[k];
    end for;
    suf[h + 1] := 1;
    for k in [2..h] do
        suf[h + 2 - k] := suf[h + 3 - k]*D[h + 2 - k];
    end for;
    Dc := Z1*suf[2];
    Xb[1] := X1*suf[2];
    for k in [2..h] do
        Xb[k] := (X1*D[k] + Delta[k])*pre[k - 1]*suf[k + 1];
    end for;
    return Xb, Dc;
end function;

function SWPushPrecomp(st, N)
    h := (N - 1) div 2;
    X1 := st[1];
    Z1 := st[2];
    Delta := st[3];
    D := st[4];
    a4 := st[5];
    a6 := st[6];
    t0 := st[7];
    // Line 7
    Xb, Dc := SWCommonDenominator(X1, Z1, Delta, D, N);
    SX := (&+[Xb[k] : k in [1..h]]);
    Dc2 := Dc^2;
    c1 := t0^2;
    return [* Xb, SX, Dc, Dc2, a4, a6, c1 *];
end function;

function SWPush(XR, ZR, pr, N)
    h := (N - 1) div 2;
    Xb := pr[1];
    SX := pr[2];
    Dc := pr[3];
    Dc2 := pr[4];
    a4 := pr[5];
    a6 := pr[6];
    c1 := pr[7];
    // Line 8
    t1 := XR*Dc;
    PiR := t1 - Xb[1]*ZR;
    K1 := 1;
    M := 0;
    for k in [2..h] do
        delta := t1 - Xb[k]*ZR;
        M := M*delta + K1;
        K1 := K1*delta + PiR;
        PiR := PiR*delta;
    end for;
    PiR2 := PiR^2;
    K2 := K1^2 - 2*PiR*M;
    t3 := XR^2;
    t5 := ZR^2;
    t4 := a4*t5;
    Fh := 4*((t3 + t4)*XR + a6*t5*ZR);
    Gh := 6*t3 + 2*t4;
    X := PiR2*(N*t1 - 2*SX*ZR) + Dc2*(Dc*Fh*K2 - Gh*K1*PiR);
    Z := ZR*Dc*PiR2;
    return c1*X, Z;
end function;

function SWCodomainOnly(a4, a6, X1, Z1, N)
    h := (N - 1) div 2;
    // Codomain-only variant of Algorithm 3 (Section 4.2): lines 1-3 without the rescaling in line 2,
    // the common denominator of line 7, and the power sums S_1, S_2, S_3.
    Fn, Gn, Delta, D, X1s, Z1s := SWMultiples(a4, a6, X1, Z1, N, 0);
    Xb, Dc := SWCommonDenominator(X1s, Z1s, Delta, D, N);
    // Power sums S_j of the kernel numerators
    S1 := (&+[Xb[k] : k in [1..h]]);
    t := AssociativeArray();
    for k in [1..h] do
        t[k] := Xb[k]^2;
    end for;
    S2 := (&+[t[k] : k in [1..h]]);
    S3 := (&+[t[k]*Xb[k] : k in [1..h]]);
    // a_4' D^4 and a_6' D^6
    t1 := Dc^2;
    t2 := t1*Dc;
    t3 := a4*t1;
    a4p := t1*((1 - 10*h)*t3 - 30*S2);
    a6p := t2*((1 - 28*h)*a6*t2 - 70*S3 - 42*S1*t3);
    return a4p, a6p;
end function;
