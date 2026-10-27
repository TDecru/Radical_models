// Sanity check of Algorithm 4 (Montgomery curves y^2 = x^3 + Ax^2 + x): three radical steps of degree N,
// compared with Magma's isogenies (codomain, next generator, pushed points), and the codomain-only variant.
// Run from this directory:  magma montgomery.m
load "lib/core.m"; load "lib/testing.m"; load "lib/montgomery.m";
SetSeed(1);
Ns := [5, 7, 9, 11, 13, 15];  // 9 and 15 are composite
p := SetupPrime(Ns, 60);
K := GF(p);
for N in Ns do
    repeat
        A := Random(K);
        E1 := EllipticCurve([K | 0, A, 0, 1, 0]);
    until #E1 mod N eq 0;
    P1 := KernelPoint(E1, N);
    An := A; C := K!1; X1 := P1[1]; Z1 := K!1;
    for step in [1..3] do
        Rs := [Random(E1) : j in [1..3]];
        Anp, Cp, Xp, Zp, st := MontgomeryStep(An, C, X1, Z1, N);
        xRs2 := [];
        for R in Rs do
            X, Z := MontgomeryPush(R[1], K!1, st, N);
            Append(~xRs2, X/Z);
        end for;
        E2 := EllipticCurve([K | 0, Anp/Cp, 0, 1, 0]);
        ok, msg := CheckStep(E1, P1, E2, Xp/Zp, N, Rs, xRs2);
        Report(Sprintf("Montgomery, N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, Xp/Zp);
        E1 := E2; An := Anp; C := Cp; X1 := Xp; Z1 := Zp;
    end for;
    Anp, Cp := MontgomeryCodomainOnly(An, C, X1, Z1, N);
    Report(Sprintf("Montgomery, N = %o, codomain only", N), CheckCodomain(E1, P1, EllipticCurve([K | 0, Anp/Cp, 0, 1, 0]), N), "");
end for;
