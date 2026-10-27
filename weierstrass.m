// Sanity check of Algorithm 3 (short Weierstrass curves y^2 = x^3 + a4 x + a6) and of the codomain-only
// variant (Velu's formulae with power sums).  Run from this directory:  magma weierstrass.m
load "lib/core.m"; load "lib/testing.m"; load "lib/weierstrass.m";
SetSeed(1);
Ns := [5, 7, 9, 11, 13, 15];  // 9 and 15 are composite
p := SetupPrime(Ns, 60);
K := GF(p);
for N in Ns do
    repeat
        a4 := Random(K); a6 := Random(K);
        E1 := EllipticCurve([K | a4, a6]);
    until #E1 mod N eq 0;
    P1 := KernelPoint(E1, N);
    X1 := P1[1];
    for step in [1..3] do
        Rs := [Random(E1) : j in [1..3]];
        a4p, a6p, xPp, st := SWStep(a4, a6, X1, K!1, N);
        pr := SWPushPrecomp(st, N);
        xRs2 := [];
        for R in Rs do
            X, Z := SWPush(R[1], K!1, pr, N);
            Append(~xRs2, X/Z);
        end for;
        E2 := EllipticCurve([K | a4p, a6p]);
        ok, msg := CheckStep(E1, P1, E2, xPp, N, Rs, xRs2);
        Report(Sprintf("short Weierstrass, N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, xPp);
        E1 := E2; a4 := a4p; a6 := a6p; X1 := xPp;
    end for;
    a4p, a6p := SWCodomainOnly(a4, a6, X1, K!1, N);
    Report(Sprintf("short Weierstrass, N = %o, codomain only", N), CheckCodomain(E1, P1, EllipticCurve([K | a4p, a6p]), N), "");
end for;
