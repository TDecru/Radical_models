// Sanity check of Algorithm 6 (curves y^2 = x(x^2 + Ax + B)) and of its codomain-only variant (trace formula).
// Run from this directory:  magma wmodel.m
load "lib/core.m"; load "lib/testing.m"; load "lib/wmodel.m";
SetSeed(1);
Ns := [5, 7, 11, 13];
p := SetupPrime(Ns, 60);
K := GF(p);
for N in Ns do
    repeat
        A := Random(K); B := Random(K);
        E1 := EllipticCurve([K | 0, A, 0, B, 0]);
    until #E1 mod N eq 0;
    P1 := KernelPoint(E1, N);
    An := A; Bn := B; C := K!1; X1 := P1[1]; Z1 := K!1;
    for step in [1..3] do
        Rs := [Random(E1) : j in [1..3]];
        Anp, Bnp, Cp, Xp, Zp, st := WStep(An, Bn, C, X1, Z1, N);
        xRs2 := [];
        for R in Rs do
            X, Z := WPush(R[1], K!1, st, N);
            Append(~xRs2, X/Z);
        end for;
        E2 := EllipticCurve([K | 0, Anp/Cp, 0, Bnp/Cp, 0]);
        ok, msg := CheckStep(E1, P1, E2, Xp/Zp, N, Rs, xRs2);
        Report(Sprintf("y^2 = x(x^2+Ax+B), N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, Xp/Zp);
        E1 := E2; An := Anp; Bn := Bnp; C := Cp; X1 := Xp; Z1 := Zp;
    end for;
    Anp, Bnp, Cp := WCodomainOnly(An, Bn, C, X1, Z1, N);
    Report(Sprintf("y^2 = x(x^2+Ax+B), N = %o, codomain only", N), CheckCodomain(E1, P1, EllipticCurve([K | 0, Anp/Cp, 0, Bnp/Cp, 0]), N), "");
end for;
