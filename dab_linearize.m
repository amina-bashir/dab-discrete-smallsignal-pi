function L = dab_linearize(D, P)
D0 = [D(1); D(2); mod(D(3),2)];
V0 = [P.V1; P.V2];
[~, G, H] = dab_period_map(zeros(5,1), D0, V0, P);
xs = (eye(5) - G)\H;
eD = 1e-6;
eV = 1e-3;
B = zeros(5,3);
for j = 1:3
    d = zeros(3,1);
    d(j) = eD;
    B(:,j) = (dab_period_map(xs, D0+d, V0, P) - dab_period_map(xs, D0-d, V0, P))/(2*eD);
end
C = zeros(5,2);
for j = 1:2
    v = zeros(2,1);
    v(j) = eV;
    C(:,j) = (dab_period_map(xs, D0, V0+v, P) - dab_period_map(xs, D0, V0-v, P))/(2*eV);
end
L.D0  = D0;
L.V0  = V0;
L.xs  = xs;
L.A   = G;
L.B   = B;
L.C   = C;
L.sys = ss(G, [B C], [0 0 1 0 0; 0 0 0 0 1], 0, P.Ts, ...
           'InputName', {'D1','D2','D3','V1','V2'}, 'OutputName', {'iLo','vCo'});
L.G   = L.sys('iLo','D3');
L.Kdc = dcgain(L.G);
end
