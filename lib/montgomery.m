// Algorithm 4 (Montgomery) and its codomain-only variant; needs core.m
// Generated from the reference implementation; names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

function MontMultiples(An, C, X1, Z1, N, withDelta)
    h := (N - 1) div 2;
    X := AssociativeArray();
    Z := AssociativeArray();
    Delta := AssociativeArray();
    // Line 1
    t1 := C*X1^2;
    t2 := C*Z1^2;
    t5 := X1*Z1;
    t3 := An*t5;
    t6 := t1 + t3 + t2;
    Fn := 0;
    if withDelta then
        Fn := 4*X1*t6;
    end if;
    Gn := 3*t1 + 2*t3 + t2;
    // Line 2
    X[1] := X1;
    Z[1] := Z1;
    X[2] := (t1 - t2)^2;
    if withDelta then
        Z[2] := C*Z1*Fn;
    else
        Z[2] := 4*C*t5*t6;
    end if;
    // Line 3
    for k in [3..h + 1] do
        if withDelta or k le h then
            t1 := (X[k - 1] - Z[k - 1])*(X1 + Z1);
            t2 := (X[k - 1] + Z[k - 1])*(X1 - Z1);
            Delta[k - 1] := (t1 - t2)/2;
            if k le h then
                X[k] := Z[k - 2]*(t1 + t2)^2;
                Z[k] := X[k - 2]*(t1 - t2)^2;
            end if;
        end if;
    end for;
    return Fn, Gn, X, Z, Delta;
end function;

function EdwardsNormalization(An, C, X, Z, N)
    // (A'_n : C') of the Costello--Hisil codomain from the twisted Edwards formula of [MR18]
    // (line 7 of Algorithm 4 if delta = 0, and the codomain-only variant).
    // Line 7
    h := (N - 1) div 2;
    t1 := (An + 2*C)^N*(&*[X[i] + Z[i] : i in [1..h]])^8;
    t2 := (An - 2*C)^N*(&*[X[i] - Z[i] : i in [1..h]])^8;
    return 2*(t1 + t2), t1 - t2;
end function;

function MontgomeryStep(An, C, X1, Z1, N)
    h := (N - 1) div 2;
    // Line 1
    Fn, Gn, X, Z, Delta := MontMultiples(An, C, X1, Z1, N, true);
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
    Pi_D := 4^(h - 2)*Z[2]*PiXpart*PiDelta[h - 1]^2;
    g, g2, g3, R, U, S := CoreRadical(Fn, C, PiDelta[h]^2, Delta, Z, N);
    Rd, FnR, abar, bbar, cbar := CoreCubic(Fn, Gn, C, Delta[2], Z[2], g, g2, g3, R, U, S);
    // Line 5
    PiX := PiXpart*X[h - 1]*X[h];
    if h eq 2 then
        PiX := X[1]*X[2];
    end if;
    q := (Z[2]*g3*Z1*Pi_D*Delta[2])^2;
    // Line 6
    t1 := PiX^2;
    t2 := t1^2;
    t3 := q^2;
    nu := t2*(9*cbar - abar*bbar) + 4*abar*t3;
    delta := t2*(6*bbar - 2*abar^2) + 12*t3;
    if delta eq 0 then
        // Line 7
        Anp, Cp := EdwardsNormalization(An, C, X, Z, N);
        Xp := t1*abar*Cp - 2*q*Anp;
        Zp := 6*q*Cp;
    else
        // Line 8
        Anp := t1*(abar*delta - 3*nu);
        Cp := 2*q*delta;
        Xp := t1*nu;
        Zp := Cp;
    end if;
    st := [* X, Z *];
    return Anp, Cp, Xp, Zp, st;
end function;

function MontgomeryPush(XR, ZR, st, N)
    h := (N - 1) div 2;
    X := st[1];
    Z := st[2];
    // Line 9
    s := AssociativeArray();
    t := AssociativeArray();
    for i in [1..h] do
        s[i] := (XR - ZR)*(X[i] + Z[i]);
        t[i] := (XR + ZR)*(X[i] - Z[i]);
    end for;
    return XR*(&*[s[i] + t[i] : i in [1..h]])^2, ZR*(&*[s[i] - t[i] : i in [1..h]])^2;
end function;

function MontgomeryCodomainOnly(An, C, X1, Z1, N)
    // Codomain-only variant of Algorithm 4 (Section 4.3): lines 1-3 restricted to X_i, Z_i,
    // and the twisted Edwards formula of line 7.
    Fn, Gn, X, Z, Delta := MontMultiples(An, C, X1, Z1, N, false);
    Anp, Cp := EdwardsNormalization(An, C, X, Z, N);
    return Anp, Cp;
end function;
