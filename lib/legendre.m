// Algorithm 10 (Legendre curves, Jacobi intersections); needs core.m, wmodel.m
// Generated from the reference implementation; names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

function LegendreStep(lamn, C, X1, Z1, N)
    h := (N - 1) div 2;
    // Line 1
    Anp, Bnp, Cp, Xp, Zp, st := WStep(-(lamn + C), lamn, C, X1, Z1, N);
    X := st[1];
    Z := st[2];
    Bn := st[3];
    Cw := st[4];
    BnN := st[5];
    CN := st[6];
    // Line 2
    pi1 := (&*[X[i] - Z[i] : i in [1..h]]);
    pil := (&*[C*X[i] - lamn*Z[i] : i in [1..h]]);
    // Line 3
    t1 := pi1^2;
    t2 := pil^2;
    lnp := BnN*CN*t1^2;
    Cp := C^2*t2^2;
    t3 := CN*t1;
    t4 := C*t2;
    // Line 4
    return lnp, Cp, t3*Xp, t4*Zp, [* st, t3, t4 *];
end function;

function LegendrePush(XR, ZR, st, N)
    stw := st[1];
    t3 := st[2];
    t4 := st[3];
    // Line 1
    X, Z := WPush(XR, ZR, stw, N);
    // Line 4
    return t3*X, t4*Z;
end function;

function LegendreCodomainOnly(lamn, C, X1, Z1, N)
    h := (N - 1) div 2;
    // Codomain-only variant of Algorithm 10: the multiples of Algorithm 6, line 2, and lambda'.
    Fn, Gn, X, Z, Delta := WMultiples(-(lamn + C), lamn, C, X1, Z1, N, false);
    // Line 2
    pi1 := (&*[X[i] - Z[i] : i in [1..h]]);
    pil := (&*[C*X[i] - lamn*Z[i] : i in [1..h]]);
    // lambda' = ((lambda_n C)^N pi_1^4 : C^2 pi_lambda^4)
    return (lamn*C)^N*pi1^4, C^2*pil^4;
end function;
