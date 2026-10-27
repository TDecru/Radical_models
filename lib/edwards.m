// Algorithm 5 (twisted Edwards); needs core.m, montgomery.m
// Generated from the reference implementation; names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

function EdwardsStep(a, d, Y1, Z1, N)
    // Line 1
    An := 2*(a + d);
    C := a - d;
    Xe := Z1 + Y1;
    Ze := Z1 - Y1;
    // Line 2
    Anp, Cp, Xp, Zp, st := MontgomeryStep(An, C, Xe, Ze, N);
    // Line 3
    return Anp + 2*Cp, Anp - 2*Cp, Xp - Zp, Xp + Zp, st;
end function;

function EdwardsPush(YR, ZR, st, N)
    // Line 1
    X, Z := MontgomeryPush(ZR + YR, ZR - YR, st, N);
    // Line 3
    return X - Z, X + Z;
end function;
