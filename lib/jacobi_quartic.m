// Algorithm 9 (extended Jacobi quartics); needs core.m, wmodel.m
// Generated from the reference implementation; names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

function JQStep(an, dn, C, W1, V1, N)
    // Line 1
    X1 := 2*(W1*C + an*V1);
    Z1 := V1*C;
    // Line 2
    Anp, Bnp, Cp, Xp, Zp, st := WStep(-4*an*C, 4*(an^2 - dn), C^2, X1, Z1, N);
    // Line 3
    anp := -Anp/4;
    dnp := (Anp^2 - 4*Bnp*Cp)/16;
    // Line 4
    return anp, dnp, Cp, Xp*Cp - 2*anp*Zp, 2*Zp*Cp, [* st, an, C, anp, Cp *];
end function;

function JQPush(WR, VR, st, N)
    stw := st[1];
    an := st[2];
    C := st[3];
    anp := st[4];
    Cp := st[5];
    // Line 1
    XR := 2*(WR*C + an*VR);
    ZR := VR*C;
    // Line 2
    X, Z := WPush(XR, ZR, stw, N);
    // Line 4
    return X*Cp - 2*anp*Z, 2*Z*Cp;
end function;

function JQCodomainOnly(an, dn, C, W1, V1, N)
    // Codomain-only variant of Algorithm 9: line 2 runs the codomain-only variant of Algorithm 6.
    // Line 1
    X1 := 2*(W1*C + an*V1);
    Z1 := V1*C;
    // Line 2
    Anp, Bnp, Cp := WCodomainOnly(-4*an*C, 4*(an^2 - dn), C^2, X1, Z1, N);
    // Line 3
    return -Anp/4, (Anp^2 - 4*Bnp*Cp)/16, Cp;
end function;
