// Algorithm 6 (y^2 = x(x^2+Ax+B)) and Proposition 4.12 (codomain only); needs core.m
// Generated from the reference implementation; names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

function WMultiples(An, Bn, C, X1, Z1, N, withDelta)
    h := (N - 1) div 2;
    X := AssociativeArray();
    Z := AssociativeArray();
    Delta := AssociativeArray();
    // Line 1
    t0 := Bn*Z1;
    X[1] := C*X1;
    Z[1] := C*Z1;
    // Line 2
    t5 := X[1]*Z[1];
    t1 := C*X[1]^2;
    t2 := Bn*Z[1]^2;
    t3 := An*t5;
    t6 := t1 + t3 + t2;
    Fn := 0;
    if withDelta then
        Fn := 4*X[1]*t6;
    end if;
    Gn := 3*t1 + 2*t3 + t2;
    X[2] := (t1 - t2)^2;
    if withDelta then
        Z[2] := C*Z[1]*Fn;
    else
        Z[2] := 4*C*t5*t6;
    end if;
    // Line 3
    for k in [3..h + 1] do
        if withDelta or k le h then
            t1 := X[k - 1]*X[1];
            t2 := Z[k - 1]*Z[1];
            Delta[k - 1] := (X[k - 1] - Z[k - 1])*(X[1] + Z[1]) - t1 + t2;
            if k le h then
                X[k] := Z[k - 2]*(t1 - t0*Z[k - 1])^2;
                Z[k] := X[k - 2]*Delta[k - 1]^2;
            end if;
        end if;
    end for;
    return Fn, Gn, X, Z, Delta;
end function;

function WTraceNormalization(An, Bn, C, X, Z, N, BnN, CN)
    // A' = A1/A2 of the codomain y^2 = x(x^2 + A'x + B^N) of phi_B (Proposition 4.12, codomain only;
    // also line 7 of Algorithm 6 if delta = 0).
    // Arithmetic in K[eta]/(eta^2 + An eta + C Bn), eta = C eps; elements are pairs [u0, u1] = u0 + u1 eta.
    // Line 7
    h := (N - 1) div 2;
    t0 := C*Bn;
    u0 := C*X[1];
    u1 := -Z[1];
    for i in [2..h] do
        v0 := C*X[i];
        v1 := -Z[i];
        w0 := u0*v0;
        w1 := u1*v1;
        u1 := (u0 + u1)*(v0 + v1) - w0 - w1 - w1*An;
        u0 := w0 - w1*t0;
    end for;
    n0 := u0^2 - u0*u1*An + u1^2*t0;
    c0 := u0 - u1*An;
    c1 := -u1;
    for j in [1..2] do
        w0 := c0^2;
        w1 := c1^2;
        c1 := 2*c0*c1 - w1*An;
        c0 := w0 - w1*t0;
    end for;
    e0 := -t0;
    e1 := -An;
    bits := Intseq(N, 2);
    if bits[#bits - 1] eq 1 then
        w0 := -e1*t0;
        e1 := e0 - e1*An;
        e0 := w0;
    end if;
    for j in [3..#bits] do
        w0 := e0^2;
        w1 := e1^2;
        e1 := 2*e0*e1 - w1*An;
        e0 := w0 - w1*t0;
        if bits[#bits + 1 - j] eq 1 then
            w0 := -e1*t0;
            e1 := e0 - e1*An;
            e0 := w0;
        end if;
    end for;
    w0 := e0*c0;
    w1 := e1*c1;
    f1 := (e0 + e1)*(c0 + c1) - w0 - w1 - w1*An;
    f0 := w0 - w1*t0;
    return -(2*f0 - f1*An), CN*n0^2;
end function;

function WStep(An, Bn, C, X1, Z1, N)
    h := (N - 1) div 2;
    // Line 1
    Fn, Gn, X, Z, Delta := WMultiples(An, Bn, C, X1, Z1, N, true);
    // Line 4
    PiXpart := X[1];
    for i in [2..h - 2] do
        PiXpart := PiXpart*X[i];
    end for;
    if h eq 2 then
        PiXpart := 1;
    end if;
    PiDelta := AssociativeArray();
    PiDelta[1] := 1;
    for k in [2..h] do
        PiDelta[k] := PiDelta[k - 1]*Delta[k];
    end for;
    Pi_D := Z[2]*PiXpart*PiDelta[h - 1]^2;
    g, g2, g3, R, U, S := CoreRadical(Fn, C, PiDelta[h]^2, Delta, Z, N);
    Rd, FnR, abar, bbar, cbar := CoreCubic(Fn, Gn, C, Delta[2], Z[2], g, g2, g3, R, U, S);
    // Line 5
    PiX := PiXpart*X[h - 1]*X[h];
    if h eq 2 then
        PiX := X[1]*X[2];
    end if;
    q := (Z[2]*g3*Z[1]*Pi_D*Delta[2])^2;
    BnN := Bn^N;
    CN := C^N;
    // Line 6
    t1 := PiX^2;
    t2 := t1^2*CN;
    t3 := q^2*BnN;
    nu := t2*(9*cbar - abar*bbar) + 4*abar*t3;
    delta := t2*(6*bbar - 2*abar^2) + 12*t3;
    if delta eq 0 then
        // Line 7
        A1, A2 := WTraceNormalization(An, Bn, C, X, Z, N, BnN, CN);
        t4 := 6*q*A2;
        Anp := 6*q*A1*CN;
        Bnp := BnN*t4;
        Cp := t4*CN;
        Xp := t1*abar*A2 - 2*q*A1;
        Zp := t4;
    else
        // Line 8
        t4 := 2*q*delta;
        Anp := t1*(abar*delta - 3*nu)*CN;
        Bnp := BnN*t4;
        Cp := t4*CN;
        Xp := t1*nu;
        Zp := t4;
    end if;
    st := [* X, Z, Bn, C, BnN, CN *];
    return Anp, Bnp, Cp, Xp, Zp, st;
end function;

function WPush(XR, ZR, st, N)
    h := (N - 1) div 2;
    X := st[1];
    Z := st[2];
    Bn := st[3];
    C := st[4];
    BnN := st[5];
    CN := st[6];
    // Line 9
    t1 := C*XR;
    t2 := t1 + Bn*ZR;
    t3 := t1 + C*ZR;
    s := AssociativeArray();
    t := AssociativeArray();
    for i in [1..h] do
        u := t1*(X[i] + Z[i]);
        s[i] := u - t2*Z[i];
        t[i] := u - t3*X[i];
    end for;
    return XR*(&*[s[i] : i in [1..h]])^2, ZR*(&*[t[i] : i in [1..h]])^2;
end function;

function WCodomainOnly(An, Bn, C, X1, Z1, N)
    // Codomain-only variant of Algorithm 6 (Proposition 4.12): lines 1-3 restricted to X_i, Z_i,
    // and the trace formula of line 7.
    Fn, Gn, X, Z, Delta := WMultiples(An, Bn, C, X1, Z1, N, false);
    BnN := Bn^N;
    CN := C^N;
    A1, A2 := WTraceNormalization(An, Bn, C, X, Z, N, BnN, CN);
    return A1*CN, BnN*A2, A2*CN;
end function;
