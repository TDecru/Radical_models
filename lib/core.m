// Algorithm 1 (radical core), shared by all models
// Generated from the reference implementation; names follow the pseudocode of the paper, and "// Line k" marks line k of the algorithm.
// Field elements live in a finite field K; curve coefficients are passed projectively (numerator, C).

// N-th root in a prime field with gcd(N, p - 1) = 1
function NthRoot(x, N)
    p := Characteristic(Parent(x));
    return x^InverseMod(N, p - 1);
end function;

function CoreRadical(Fn, C, PD2, Delta, D, N)
    // Lines 1-8 of Algorithm 1. Algorithm 1 is split into two functions: this one (lines 1-8) and
    // CoreCubic (lines 9-10). The input G_n of Algorithm 1 only enters at line 9, so it is an input of
    // CoreCubic, not of this function. Algorithm 2 (Tate normal form) only calls this function; all
    // other models call both, in this order.
    // PD2 = Pi_Delta^2, supplied by the caller (Lemma 3.6). The radicand uses
    // Pi_D^2 mu^(N-1) = (mu Pi_eta[h])^2, so neither Pi_D nor mu^(N-1) is needed.
    h := (N - 1) div 2;
    // Line 1
    mu := Fn*PD2;
    eta := AssociativeArray();
    for k in [2..h] do
        eta[k] := mu*D[k];
    end for;
    // Line 2
    Pi_eta := AssociativeArray();
    Pi_eta[1] := 1;
    for i in [2..h] do
        Pi_eta[i] := Pi_eta[i - 1]*eta[i];
    end for;
    gamma := NthRoot(-C*(mu*Pi_eta[h])^2, N);
    // Line 3
    xi := AssociativeArray();
    for k in [2..h] do
        xi[k] := gamma*Delta[k];
    end for;
    for k in [h + 1..N - 2] do
        xi[k] := xi[N - k];
        eta[k] := eta[N - k];
    end for;
    // Line 4
    // Running products Pi_xi[i] = gamma xi_2 ... xi_i (they include the extra factor gamma) and Pi_eta[j] = eta_2 ... eta_j; e_i = Pi_xi[i] Pi_eta[N-1-i]
    Pi_xi := AssociativeArray();
    e := AssociativeArray();
    Pi_xi[1] := gamma;
    for i in [2..N - 2] do
        Pi_xi[i] := Pi_xi[i - 1]*xi[i];
    end for;
    for i in [h + 1..N - 2] do
        Pi_eta[i] := Pi_eta[i - 1]*eta[i];
    end for;
    for i in [1..N - 2] do
        e[i] := Pi_xi[i]*Pi_eta[N - 1 - i];
    end for;
    // Line 5
    Lam := e[1]*e[N - 2];
    p := AssociativeArray();
    for i in [1..N - 3] do
        p[i] := e[i]*e[i + 1];
    end for;
    s := AssociativeArray();
    u := AssociativeArray();
    w := AssociativeArray();
    for i in [1..h - 1] do
        s[i] := e[i] + e[N - 1 - i];
        u[i] := e[i] + e[N - 2 - i];
    end for;
    for i in [1..h - 2] do
        w[i] := e[i] + e[N - 3 - i];
    end for;
    Sigma1 := (&+[s[i] : i in [1..h - 1]]) + e[h];
    Sigma11 := (&+[s[i]^2 : i in [1..h - 1]]) - (N - 4)*Lam;
    Sigma2 := (&+[p[i] : i in [1..N - 3]]);
    if h gt 2 then
        Sigma3 := (&+[w[i]*(p[i + 1] + p[N - 2 - i]) : i in [1..h - 2]]) - Lam*((&+[s[i] : i in [2..h - 1]]) - e[h]);
    else
        Sigma3 := Lam*e[h];
    end if;
    Sigma12 := (&+[u[i]*(p[i] + p[N - 2 - i]) : i in [1..h - 1]]) - Lam*(&+[u[i] : i in [1..h - 1]]);
    // Line 6
    g := mu*e[1];
    // Line 7
    g2 := g^2;
    t1 := mu*Sigma2;
    V := mu^2*Sigma3 - mu*g*e[N - 2]^2 - g2*e[N - 2];
    T1 := g*(Sigma1 - e[N - 2]) + t1;
    g3 := g2*g;
    T2 := g3 - mu*V;
    // Line 8
    t2 := D[2]*T2;
    U := t2 + 3*xi[2]*Sigma1*g2;
    S := xi[2]*(g*(2*Sigma1 + e[N - 2]) - t1);
    R := 2*xi[2]*T1*t2 + xi[2]^2*g2*(3*(Sigma1^2 - Sigma11)*g - 2*mu*Sigma12);
    return g, g2, g3, R, U, S;
end function;

function CoreCubic(Fn, Gn, C, Delta2, D2, g, g2, g3, R, U, S)
    // Lines 9-10 of Algorithm 1, continuing CoreRadical (lines 1-8): W, R_d, B and the scaled cubic
    // (abar, bbar, cbar) of Proposition 3.5. This is where the input G_n of Algorithm 1 is used
    // (in W = F_n S - G_n Delta_2 g^2 and in B). Delta2 and D2 are Delta_2 and D_2, and g, g2, g3,
    // R, U, S are the outputs of CoreRadical.
    // Line 9
    t1 := Gn*Delta2;
    t2 := Fn*S;
    W := t2 - t1*g2;
    Rd := 2*Delta2^2*g3*g2;
    FnR := Fn*R;
    B := C*Delta2*Rd + D2*(2*t2*S*g - 4*t1*S*g3 - 3*FnR);
    // Line 10
    t3 := D2*g;
    abar := g*(2*t3*W^2 - Fn*B);
    t4 := B*U;
    t5 := Fn^2;
    t6 := t3*t5*g;
    bbar := -2*t6*t4*W;
    cbar := t6*t5*t4^2/2;
    return Rd, FnR, abar, bbar, cbar;
end function;
