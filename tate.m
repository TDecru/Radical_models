// Sanity check of Algorithm 2 (Tate normal form y^2 + (1-c)xy - by = x^3 - bx^2, P = (0,0)).
// Run from this directory:  magma tate.m
load "lib/core.m"; load "lib/testing.m"; load "lib/tate.m";
SetSeed(1);
Ns := [5, 7, 11, 13];
p := SetupPrime(Ns, 60);
K := GF(p);
// (b, c) of (E, P) for E: y^2 = x^3 + a4 x + a6 (Lemma 3.2)
function TateNormalForm(E, P)
    a := Coefficients(E);
    x1 := P[1];
    F := 4*(x1^3 + a[4]*x1 + a[5]);
    G := 6*x1^2 + 2*a[4];
    xh2 := (2*P)[1] - x1;
    return xh2^3/F, 1 + G*xh2/F;
end function;
for N in Ns do
    repeat
        E0 := EllipticCurve([K | Random(K), Random(K)]);
    until #E0 mod N eq 0;
    b, c := TateNormalForm(E0, KernelPoint(E0, N));
    bn := b; cn := c; C := K!1;
    for step in [1..3] do
        E1 := EllipticCurve([K | 1 - cn/C, -bn/C, -bn/C, 0, 0]);
        P1 := E1![0, 0];
        Rs := [Random(E1) : j in [1..3]];
        bnp, cnp, Cp, st := TateStep(bn, cn, C, N);
        pr := TatePushPrecomp(st, N);
        xRs2 := [];
        for R in Rs do
            X, Z := TatePush(C^2*R[1], K!1, pr, N);     // the algorithm works in the model scaled by C
            Append(~xRs2, X/(Z*Cp^2));
        end for;
        E2 := EllipticCurve([K | 1 - cnp/Cp, -bnp/Cp, -bnp/Cp, 0, 0]);
        ok, msg := CheckStep(E1, P1, E2, K!0, N, Rs, xRs2);
        Report(Sprintf("Tate normal form, N = %o, step %o", N, step), ok, msg);
        bn := bnp; cn := cnp; C := Cp;
    end for;
end for;
