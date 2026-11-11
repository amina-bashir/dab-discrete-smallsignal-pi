clear; clc;
[P, M] = dab_params();
k = 1;

D1 = M.D(k,1);
D2 = M.D(k,2);
D3 = M.D(k,3);
Vi = P.V1;  Vo = P.V2;
Ri = P.R1;  Ro = P.R2;
Li = P.L1;  Lo = P.L2;
Ci = P.C1;  Co = P.C2;
Ls = P.L;   Rs = P.R;
n  = P.n;   f_sw = P.fs;
Rds = 1e-6;  Vth = 0.5;  G_mosfet = 1e-9;
Tsolver = 0.5e-6;
tdead = 0;

kp_cmd = -5.2995e-05;
ki_cmd = -5.2995;
umax_cmd   = 0.8*min(D3 - 0, (D1 - D2) - D3);
i2_cmd_max =  umax_cmd;
i2_cmd_min = -umax_cmd;

Iref0 = -41.3;
Iref1 = -46.3;
t_step = 5e-3;

fprintf('%s  D = [%.2f %.2f %.2f]\n', M.label{k}, D1, D2, D3);
fprintf('Kp = %.4g   Ki = %.4g   clamp on D3 = +/- %.3f\n', kp_cmd, ki_cmd, umax_cmd);
fprintf('Reference step %.1f A -> %.1f A at %.0f ms\n', Iref0, Iref1, t_step*1e3);
