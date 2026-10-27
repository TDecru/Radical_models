// Sanity check of Algorithm 7 (Huff curves a x(y^2 - 1) = b y(x^2 - 1), coordinate m = y/x).
// The comparison uses V^2 = U(U + a^2)(U + b^2) with U = ab(am - b)/(a - bm) (Joye, Tibouchi and Vergnaud).
// Run from this directory:  magma huff.m
load "lib/core.m"; load "lib/testing.m"; load "lib/montgomery.m"; load "lib/huff.m";
SetSeed(1);
Ns := [5, 7, 9, 11, 13, 15];  // 9 and 15 are composite
p := SetupPrime(Ns, 60);
K := GF(p);
function HuffW(a, b)
    return EllipticCurve([Parent(a) | 0, a^2 + b^2, 0, a^2*b^2, 0]);
end function;
function HuffU(a, b, M1, M2)
    m := M1/M2;
    return a*b*(a*m - b)/(a - b*m);
end function;
for N in Ns do
    repeat
        a := Random(K); b := Random(K);
        E1 := HuffW(a, b);
    until #E1 mod N eq 0;
    P1 := KernelPoint(E1, N);
    U := P1[1]; M1 := a*(U + b^2); M2 := b*(U + a^2);
    for step in [1..3] do
        Rs := [Random(E1) : j in [1..3]];
        ap, bp, M1p, M2p, st := HuffStep(a, b, M1, M2, N);
        xRs2 := [];
        for R in Rs do
            U := R[1];
            Mr1, Mr2 := HuffPush(a*(U + b^2), b*(U + a^2), st, N);
            Append(~xRs2, HuffU(ap, bp, Mr1, Mr2));
        end for;
        E2 := HuffW(ap, bp);
        ok, msg := CheckStep(E1, P1, E2, HuffU(ap, bp, M1p, M2p), N, Rs, xRs2);
        Report(Sprintf("Huff, N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, HuffU(ap, bp, M1p, M2p));
        E1 := E2; a := ap; b := bp; M1 := M1p; M2 := M2p;
    end for;
    ap, bp := HuffCodomainOnly(a, b, M1, M2, N);
    Report(Sprintf("Huff, N = %o, codomain only", N), CheckCodomain(E1, P1, HuffW(ap, bp), N), "");
end for;
