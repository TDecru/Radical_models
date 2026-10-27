// Algorithm 11 (twisted Hessian curves); needs core.m
// Names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

function HessMultiples(a, d, Th, Xi, N, rs)
    // Lines 1-4. If rs = 1, line 3 rescales (X_1 : Z_1) by F_n, so that D_2 = 1.
    h := (N - 1) div 2;
    Delta := AssociativeArray();
    D := AssociativeArray();
    // Line 1
    d2 := d^2;
    k := d2*d - 27*a;
    X1 := k*Th;
    Z1 := d*Th + 3*Xi;
    // Line 2
    t1 := X1^2;
    t2 := Z1^2;
    t3 := d2*X1*Z1;
    kd := k*d;
    t4 := kd*t2;
    t5 := k*Z1;
    t6 := t5^2*Z1;
    Fn := X1*(4*t1 - 9*t3 + 6*t4) - t6;
    Gn := (6*t1 - 9*t3 + 3*t4)/2;
    // Line 3
    Delta[1] := 0;
    D[1] := 1;
    Delta[2] := 3*(X1*t6 - t1*(t1 - 3*t3 + 3*t4));
    D[2] := Fn;
    if rs eq 1 then
        t7 := Fn^2;
        X1 := X1*Fn;
        Z1 := Z1*Fn;
        t1 := t1*t7;
        t5 := t5*Fn;
        Gn := Gn*t7;
        Fn := t7^2;
        D[2] := 1;
    end if;
    // Line 4
    for j in [2..h - 1] do
        t7 := Delta[j]^2;
        Delta[j + 1] := D[j - 1]*D[j]*(Fn*D[j] + 2*Gn*Delta[j]) - t7*Delta[j - 1];
        D[j + 1] := t7*D[j - 1];
    end for;
    return d2, k, kd, X1, Z1, t1, t5, Fn, Gn, Delta, D;
end function;

function HessKernelSums(a, d, k, X1, t1, t5, Delta, D, N)
    // Pi = c^h Pi_n / Pi_m and Sigma = c Sigma_n / Pi_m (codomain only, and the exceptional case).
    h := (N - 1) div 2;
    Xt := AssociativeArray();
    Xt2 := AssociativeArray();
    n := AssociativeArray();
    m := AssociativeArray();
    // Line 8
    Xt[1] := X1;
    Xt2[1] := t1;
    n[1] := t1;
    m[1] := (t5 - d*X1)^3 + 27*a*X1*t1;
    Pi_m := m[1];
    Pi_n := n[1];
    Sigma_n := n[1];
    for i in [2..h] do
        Xt[i] := X1*D[i] + Delta[i];
        t7 := Xt[i]^2;
        Xt2[i] := t7;
        n[i] := t7*D[i];
        m[i] := (t5*D[i] - d*Xt[i])^3 + 27*a*t7*Xt[i];
        Sigma_n := Sigma_n*m[i] + Pi_m*n[i];
        Pi_m := Pi_m*m[i];
        Pi_n := Pi_n*n[i];
    end for;
    // Line 8
    c := 27*t5;
    dn := N*d*Pi_m - 6*a*c*Sigma_n;
    dd := c^h*Pi_n;
    return Xt, Xt2, Pi_m, Pi_n, dn, dd;
end function;

function HessKernelProduct(a, d, X1, t1, t5, Delta, D, Pi_D, N)
    // Line 6: Pi = c^h Pi_n / Pi_m, with Pi_n = (prod Xt_i)^2 Pi_D.
    h := (N - 1) div 2;
    Xt := AssociativeArray();
    Xt2 := AssociativeArray();
    m := AssociativeArray();
    // Line 6
    Xt[1] := X1;
    Xt2[1] := t1;
    m[1] := (t5 - d*X1)^3 + 27*a*X1*t1;
    Pi_m := m[1];
    PX := X1;
    for i in [2..h] do
        Xt[i] := X1*D[i] + Delta[i];
        t7 := Xt[i]^2;
        m[i] := (t5*D[i] - d*Xt[i])^3 + 27*a*t7*Xt[i];
        Xt2[i] := t7;
        PX := PX*Xt[i];
        Pi_m := Pi_m*m[i];
    end for;
    Pi_n := PX^2*Pi_D;
    return Xt, Xt2, Pi_m, Pi_n;
end function;

function HessianStep(a, d, Th, Xi, N)
    h := (N - 1) div 2;
    // Line 1
    d2, k, kd, X1, Z1, t1, t5, Fn, Gn, Delta, D := HessMultiples(a, d, Th, Xi, N, 1);
    // Line 5
    Pi_D := (&*[D[k] : k in [2..h]]);
    PD2 := Delta[h]^2*D[h - 1]*D[h];
    g, g2, g3, R, U, S := CoreRadical(Fn, 1, PD2, Delta, D, N);
    Rd, FnR, abar, bbar, cbar := CoreCubic(Fn, Gn, 1, Delta[2], D[2], g, g2, g3, R, U, S);
    // Line 6
    Xt, Xt2, Pi_m, Pi_n := HessKernelProduct(a, d, X1, t1, t5, Delta, D, Pi_D, N);
    // Line 7
    c := 27*t5;
    dd := c^h*Pi_n;
    T := Rd*g*Z1*Fn;
    tau := dd*T;
    t2 := abar^2;
    A4 := t2 - 3*bbar;
    A6 := abar*(2*t2 - 9*bbar) + 27*cbar;
    t3 := Pi_m^2;
    t4 := t3^2;
    t6 := T^2;
    C4 := 16*A4*t4*t6/9;
    C6 := -32*A6*t4*t3*t6*T/27;
    aN := a^N;
    atau := aN*tau^3;
    t7 := atau^2;
    C6p := 157464*t7 + C6;
    C6m := 577368*t7 - C6;
    t9 := C4^2*C4;
    t10 := atau*C6p;
    nu := 756*atau*C4*(C6p*C6m + t9);
    eps := 571536*t10^2 + t9*(C6p + C6m);
    // Line 8
    if eps eq 0 then
        Xt, Xt2, Pi_m, Pi_n, dn, dd := HessKernelSums(a, d, k, X1, t1, t5, Delta, D, N);
        ap := aN*dd^3;
        dp := dn;
        t8 := Pi_m^2;
    else
        ap := atau*eps^3;
        dp := nu;
        t8 := t3*(T*eps)^2;
    end if;
    // Line 9
    XV := X1*Rd + FnR;
    ZV := Z1*Rd;
    // Line 10
    t9 := dp^2;
    t10 := t8*(4*XV - 3*d2*ZV) + 3*t9*ZV;
    t11 := 4*ZV;
    // Line 11
    kp := t9*dp - 27*ap;
    st := [* d, k, kd, Z1, t5, t8, Pi_D, Xt, Xt2, D, Pi_n, kp, dp *];
    return ap, dp, 3*t10, kp*t11 - dp*t10, st;
end function;

function HessianPushPrecomp(st, N)
    h := (N - 1) div 2;
    d := st[1];
    k := st[2];
    kd := st[3];
    Z1 := st[4];
    t5 := st[5];
    t8 := st[6];
    Pi_D := st[7];
    Xt := st[8];
    Xt2 := st[9];
    D := st[10];
    Pi_n := st[11];
    kp := st[12];
    dp := st[13];
    // Line 12
    v1 := 3*kd*Z1;
    v2 := -t5*k;
    v3 := t8*Z1^(N - 1)*Pi_D^3;
    v4 := v2*Z1;
    XtD := AssociativeArray();
    D2 := AssociativeArray();
    for i in [1..h] do
        XtD[i] := Xt[i]*D[i];
        D2[i] := D[i]^2;
    end for;
    return [* v1, v2, v3, v4, XtD, D2 *];
end function;

function HessianPush(ThR, XiR, st, pr, N)
    h := (N - 1) div 2;
    d := st[1];
    k := st[2];
    kd := st[3];
    Z1 := st[4];
    t5 := st[5];
    t8 := st[6];
    Pi_D := st[7];
    Xt := st[8];
    Xt2 := st[9];
    D := st[10];
    Pi_n := st[11];
    kp := st[12];
    dp := st[13];
    v1 := pr[1];
    v2 := pr[2];
    v3 := pr[3];
    v4 := pr[4];
    XtD := pr[5];
    D2 := pr[6];
    // Line 13
    XR := k*ThR;
    ZR := d*ThR + 3*XiR;
    w0 := XR^2;
    w1 := XR*ZR;
    w2 := -(v1*w1 + v2*ZR^2);
    w3 := -v4*w1;
    w4 := XR*Z1;
    X := v3*XR*(&*[w0*Xt2[i] + w2*XtD[i] + w3*D2[i] : i in [1..h]]);
    Z := ZR*Pi_n*(&*[w4*D[i] - ZR*Xt[i] : i in [1..h]])^2;
    return 3*X, kp*Z - dp*X;
end function;

function HessianCodomainOnly(a, d, Th, Xi, N)
    // Codomain-only variant of Algorithm 11: lines 1-4 without the rescaling in line 3,
    // and Theorem 4.26 with the sums of line 8.
    // Line 1
    d2, k, kd, X1, Z1, t1, t5, Fn, Gn, Delta, D := HessMultiples(a, d, Th, Xi, N, 0);
    // Line 8
    Xt, Xt2, Pi_m, Pi_n, dn, dd := HessKernelSums(a, d, k, X1, t1, t5, Delta, D, N);
    // a' = a^N d_d^3 and d' = d_n (Theorem 4.26)
    return a^N*dd^3, dn;
end function;
