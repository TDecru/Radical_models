// Algorithm 8 (Kummer lines); needs core.m, montgomery.m, huff.m
// Names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

function KummerStep(a2, b2, X, Z, N)
    // Line 1
    a2p, b2p, Xp, Zp, st := HuffStep(a2, b2, X, Z, N);
    // Line 2
    return a2p, b2p, Xp, Zp, st;
end function;

function KummerPush(XR, ZR, st, N)
    // Line 1
    X, Z := HuffPush(XR, ZR, st, N);
    // Line 2
    return X, Z;
end function;
