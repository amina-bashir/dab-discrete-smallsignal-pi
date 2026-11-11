function [P, M] = dab_params()
P.V1 = 800;
P.V2 = 800;
P.R1 = 100e-3;
P.R2 = 100e-3;
P.L1 = 10e-6;
P.L2 = 10e-6;
P.C1 = 30e-6;
P.C2 = 44e-6;
P.L  = 15e-6;
P.R  = 25e-3;
P.n  = 1.7;
P.fs = 50e3;
P.Ts = 1/P.fs;
P.w  = 2*pi*P.fs;

M.name  = {'fwd1','fwd2','fwd3','fwd4','fwd5','fwd6', ...
           'rev1','rev2','rev3','rev4','rev5','rev6'};
M.label = {'Mode 1 G2V','Mode 2 G2V','Mode 3 G2V','Mode 4 G2V','Mode 5 G2V','Mode 6 G2V', ...
           'Mode 1'' V2G','Mode 2'' V2G','Mode 3'' V2G','Mode 4'' V2G','Mode 5'' V2G','Mode 6'' V2G'};
M.D = [0.80 0.30  0.40
       0.30 0.60  0.80
       0.35 0.30  0.50
       0.40 0.50  0.60
       0.70 0.40  0.35
       0.80 0.70  0.50
       0.80 0.30 -0.70
       0.60 0.80 -0.10
       0.35 0.30 -0.40
       0.40 0.50 -0.45
       0.70 0.40 -0.50
       0.80 0.70 -0.40];
end
