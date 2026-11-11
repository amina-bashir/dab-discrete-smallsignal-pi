function [x, G, H, t, X] = dab_period_map(x, D, Vg, P, ns)
if nargin < 5
    ns = 0;
end
[tk, s1, s2] = tps_edges(D, P.w);
f = [0; Vg(1)/P.L1; Vg(2)/P.L2; 0; 0];
G = eye(5);
H = zeros(5,1);
t = 0;
X = x;
for k = 1:8
    A = dab_A(s1(k), s2(k), P);
    if ns > 0
        [Ed, pd] = interval_map(A, f, tk(k)/ns);
        xi = X(:,end);
        for q = 1:ns
            xi = Ed*xi + pd;
            t(end+1) = t(end) + tk(k)/ns;
            X(:,end+1) = xi;
        end
    end
    [E, p] = interval_map(A, f, tk(k));
    x = E*x + p;
    G = E*G;
    H = E*H + p;
end
end

function [tk, s1, s2] = tps_edges(D, w)
ph = pi*D(:)';
e = sort([0, mod(ph(1),2*pi), pi, mod(pi+ph(1),2*pi), mod(ph(3),2*pi), ...
          mod(ph(3)+ph(2),2*pi), mod(pi+ph(3),2*pi), mod(pi+ph(3)+ph(2),2*pi)]);
tk = zeros(1,8);
s1 = tk;
s2 = tk;
for k = 1:8
    lo = e(k);
    hi = e(mod(k,8)+1);
    if k == 8
        hi = hi + 2*pi;
    end
    m = (lo + hi)/2;
    s1(k) = 0.5*(sign(sin(m)) - sign(sin(m - ph(1))));
    s2(k) = 0.5*(sign(sin(m - ph(3))) - sign(sin(m - ph(3) - ph(2))));
    tk(k) = (hi - lo)/w;
end
end

function A = dab_A(s1, s2, P)
A = [ -P.R/P.L,       0,           0,           s1/P.L,  -P.n*s2/P.L
       0,            -P.R1/P.L1,   0,          -1/P.L1,   0
       0,             0,          -P.R2/P.L2,   0,       -1/P.L2
      -s1/P.C1,       1/P.C1,      0,           0,        0
       P.n*s2/P.C2,   0,           1/P.C2,      0,        0 ];
end

function [E, p] = interval_map(A, f, dt)
M = expm([A, f; zeros(1,6)]*dt);
E = M(1:5,1:5);
p = M(1:5,6);
end
