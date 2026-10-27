// Algorithm 2 (Tate normal form); needs core.m
// Generated from the reference implementation; names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

function TateStep(bn, cn, C, N)
    h := (N - 1) div 2;
    // Line 1
    X := AssociativeArray();
    D := AssociativeArray();
    a1 := C - cn;
    X[2] := bn*C;
    X[3] := cn*C;
    a3 := -X[2]*C;
    D[2] := One(Parent(C));
    D[3] := One(Parent(C));
    // Line 2
    for k in [4..h] do
        t1 := a3*D[k - 1];
        t2 := X[k - 1]^2;
        X[k] := t1*(t1 + a1*X[k - 1])*D[k - 2] - X[k - 2]*t2;
        D[k] := t2*D[k - 2];
    end for;
    // Line 3
    PD2 := X[2]^2;
    if h ge 3 then
        PD2 := PD2*X[h]^2*D[h - 1]*D[h];
    end if;
    g, g2, g3, R, U, S := CoreRadical(a3^2, 1, PD2, X, D, N);
    // Line 4
    t1 := C*S;
    t2 := C*U;
    bnp := g*(2*bn*g3*g2 + 2*a1*S*g3 + 2*t1*S*g - 3*C*R);
    cnp := 2*U*(t2 - g*(a1*g2 + 2*t1));
    Cp := 2*t2*U;
    st := [* X, D, a1, a3, bn, C, g, g2, g3, R, U *];
    return bnp, cnp, Cp, st;
end function;

function TatePushPrecomp(st, N)
    h := (N - 1) div 2;
    X := st[1];
    D := st[2];
    a1 := st[3];
    a3 := st[4];
    bn := st[5];
    C := st[6];
    g := st[7];
    g2 := st[8];
    g3 := st[9];
    R := st[10];
    U := st[11];
    // Line 5
    pre := AssociativeArray();
    suf := AssociativeArray();
    Xb := AssociativeArray();
    pre[3] := 1;
    for k in [4..h - 1] do
        pre[k] := pre[k - 1]*D[k];
    end for;
    suf[4] := 1;
    suf[h + 1] := 1;
    for k in [4..h] do
        suf[h + 4 - k] := suf[h + 5 - k]*D[h + 4 - k];
    end for;
    Dc := suf[4];
    Xb[1] := Zero(Parent(C));
    Xb[2] := X[2]*Dc;
    Xb[3] := X[3]*Dc;
    for k in [4..h] do
        Xb[k] := X[k]*pre[k - 1]*suf[k + 1];
    end for;
    SX := (&+[Xb[k] : k in [1..h]]);
    Dc2 := Dc^2;
    b2 := a1^2 - 4*bn*C;
    b4 := a1*a3;
    b6 := a3^2;
    v1 := 2*g*U^2;
    v2 := 2*g3*g2;
    v3 := C^2*R;
    return [* Xb, SX, Dc, Dc2, b2, b4, b6, v1, v2, v3 *];
end function;

function TatePush(XR, ZR, pr, N)
    h := (N - 1) div 2;
    Xb := pr[1];
    SX := pr[2];
    Dc := pr[3];
    Dc2 := pr[4];
    b2 := pr[5];
    b4 := pr[6];
    b6 := pr[7];
    v1 := pr[8];
    v2 := pr[9];
    v3 := pr[10];
    // Line 6
    t1 := XR*Dc;
    PiR := t1;
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
    t3 := ZR^2;
    t4 := b4*t3;
    t5 := (4*XR + b2*ZR)*XR;
    Fh := (t5 + 2*t4)*XR + b6*t3*ZR;
    Gh := t5 + 2*XR^2 + t4;
    X1 := PiR2*(N*t1 - 2*SX*ZR) + Dc2*(Dc*Fh*K2 - Gh*K1*PiR);
    Z1 := ZR*Dc*PiR2;
    return v1*(v2*X1 - v3*Z1), Z1;
end function;
