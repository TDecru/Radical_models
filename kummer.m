// Sanity check of Algorithm 8 (Kummer lines K_{a^2,b^2} in squared theta coordinates [X : Z]).
// The comparison with Velu uses the Montgomery curve with A = a^2/b^2 + b^2/a^2 and
// x = (a^2 X - b^2 Z)/(a^2 Z - b^2 X) (Proposition 4.16).  Independently of that model, the Kummer line arithmetic
// of Karati and Sarkar (Table 2 of their paper) checks that the conversion respects doubling and differential addition, that the kernel
// generator and the next generator have order N on their Kummer lines, and that pushing commutes with doubling.
// Run from this directory:  magma kummer.m
load "lib/core.m"; load "lib/testing.m"; load "lib/montgomery.m"; load "lib/huff.m"; load "lib/kummer.m";
SetSeed(1);
Ns := [5, 7, 9, 11, 13, 15];  // 9 and 15 are composite
p := SetupPrime(Ns, 60);
K := GF(p);
function KummerW(a2, b2)
    return EllipticCurve([Parent(a2) | 0, a2/b2 + b2/a2, 0, 1, 0]);
end function;
function KummerX(a2, b2, X, Z)
    return (a2*X - b2*Z)/(a2*Z - b2*X);
end function;
// Kummer line arithmetic of Karati and Sarkar (Table 2) on K_{a2,b2}, with identity [a2 : b2]; points are pairs [X, Z].
function KDbl(a2, b2, P)
    s0 := (a2 - b2)*(P[1] + P[2])^2; t0 := (a2 + b2)*(P[1] - P[2])^2;
    return [b2*(s0 + t0)^2, a2*(s0 - t0)^2];
end function;
function KAdd(a2, b2, P, Q, D)                  // P + Q, given D = P - Q
    s0 := (a2 - b2)*(P[1] + P[2])*(Q[1] + Q[2]); t0 := (a2 + b2)*(P[1] - P[2])*(Q[1] - Q[2]);
    return [D[2]*(s0 + t0)^2, D[1]*(s0 - t0)^2];
end function;
function KMul(a2, b2, n, P)                     // [n]P by the Montgomery ladder, n >= 1
    R0 := P; R1 := KDbl(a2, b2, P);
    bits := Reverse(Intseq(n, 2));
    for i in [2..#bits] do
        if bits[i] eq 1 then R0 := KAdd(a2, b2, R1, R0, P); R1 := KDbl(a2, b2, R1);
        else R1 := KAdd(a2, b2, R1, R0, P); R0 := KDbl(a2, b2, R0); end if;
    end for;
    return R0;
end function;
function KEq(P, Q)
    return P[1]*Q[2] eq P[2]*Q[1];
end function;
// [X : Z] has order exactly N on K_{a2,b2}.
function KOrder(a2, b2, P, N)
    I := [a2, b2];
    return KEq(KMul(a2, b2, N, P), I) and &and[not KEq(KMul(a2, b2, N div q, P), I) : q in PrimeDivisors(N)];
end function;
for N in Ns do
    repeat
        a2 := Random(K); b2 := Random(K);
        E1 := KummerW(a2, b2);
    until #E1 mod N eq 0;
    P1 := KernelPoint(E1, N);
    X1 := a2*P1[1] + b2; Z1 := b2*P1[1] + a2;
    for step in [1..3] do
        Rs := [Random(E1) : j in [1..3]];
        PK := [X1, Z1];
        okK := KummerX(a2, b2, KDbl(a2, b2, PK)[1], KDbl(a2, b2, PK)[2]) eq (2*P1)[1];
        PK3 := KAdd(a2, b2, KDbl(a2, b2, PK), PK, PK);
        okK := okK and KummerX(a2, b2, PK3[1], PK3[2]) eq (3*P1)[1] and KOrder(a2, b2, PK, N);
        a2p, b2p, Xp, Zp, st := KummerStep(a2, b2, X1, Z1, N);
        okK := okK and KOrder(a2p, b2p, [Xp, Zp], N);
        xRs2 := [];
        for R in Rs do
            X, Z := KummerPush(a2*R[1] + b2, b2*R[1] + a2, st, N);
            Append(~xRs2, KummerX(a2p, b2p, X, Z));
            if (2*N)*R ne E1!0 then
                RK2 := KDbl(a2, b2, [a2*R[1] + b2, b2*R[1] + a2]);
                X2, Z2 := KummerPush(RK2[1], RK2[2], st, N);
                okK := okK and KEq(KDbl(a2p, b2p, [X, Z]), [X2, Z2]);
            end if;
        end for;
        E2 := KummerW(a2p, b2p);
        ok, msg := CheckStep(E1, P1, E2, KummerX(a2p, b2p, Xp, Zp), N, Rs, xRs2);
        if not okK then ok := false; msg := "Kummer line arithmetic (Karati--Sarkar)"; end if;
        Report(Sprintf("Kummer line, N = %o, step %o", N, step), ok, msg);
        _, P1 := IsPoint(E2, KummerX(a2p, b2p, Xp, Zp));
        E1 := E2; a2 := a2p; b2 := b2p; X1 := Xp; Z1 := Zp;
    end for;
end for;
