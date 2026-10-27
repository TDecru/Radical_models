// Sanity check of Algorithm 9 (extended Jacobi quartics y^2 = d x^4 + 2a x^2 + 1, coordinate w = (1+y)/x^2).
// The comparison uses the model v^2 = u(u^2 - 4a u + 4(a^2 - d)) with u = 2w + 2a (Lemma 4.18).
// Run from this directory:  magma jacobi_quartic.m
load "lib/core.m"; load "lib/testing.m"; load "lib/wmodel.m"; load "lib/jacobi_quartic.m";
SetSeed(1);
Ns := [5, 7, 9, 11, 13, 15];  // 9 and 15 are composite
p := SetupPrime(Ns, 60);
K := GF(p);
function JQW(a, d)
    return EllipticCurve([Parent(a) | 0, -4*a, 0, 4*(a^2 - d), 0]);
end function;
for N in Ns do
    repeat
        a := Random(K); d := Random(K);
        E1 := JQW(a, d);
    until #E1 mod N eq 0;
    P1 := KernelPoint(E1, N);
    an := a; dn := d; C := K!1; W1 := (P1[1] - 2*a)/2; V1 := K!1;
    for step in [1..3] do
        Rs := [Random(E1) : j in [1..3]];
        anp, dnp, Cp, Wp, Vp, st := JQStep(an, dn, C, W1, V1, N);
        ap := anp/Cp; dp := dnp/Cp^2;
        xRs2 := [];
        for R in Rs do
            Wr, Vr := JQPush((R[1] - 2*an/C)/2, K!1, st, N);
            Append(~xRs2, 2*Wr/Vr + 2*ap);
        end for;
        E2 := JQW(ap, dp);
        ok, msg := CheckStep(E1, P1, E2, 2*Wp/Vp + 2*ap, N, Rs, xRs2);
        Report(Sprintf("Jacobi quartic, N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, 2*Wp/Vp + 2*ap);
        E1 := E2; an := anp; dn := dnp; C := Cp; W1 := Wp; V1 := Vp;
    end for;
    anp, dnp, Cp := JQCodomainOnly(an, dn, C, W1, V1, N);
    Report(Sprintf("Jacobi quartic, N = %o, codomain only", N), CheckCodomain(E1, P1, JQW(anp/Cp, dnp/Cp^2), N), "");
end for;
