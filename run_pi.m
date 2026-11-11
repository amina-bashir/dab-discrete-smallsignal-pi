clear; close all; clc;
root = fileparts(mfilename('fullpath'));
rdir = fullfile(root, 'results');
if ~isfolder(rdir), mkdir(rdir); end

[P, M] = dab_params();
nm = numel(M.name);
fc = 300;
N  = 1;
opt = pidtuneOptions('PhaseMargin', 60);
zd = tf(1, [1 zeros(1,N)], P.Ts);

G   = cell(nm,1);
Kdc = zeros(nm,1);
Kpn = zeros(nm,1);
Kin = zeros(nm,1);
for k = 1:nm
    L = dab_linearize(M.D(k,:), P);
    G{k} = L.G*zd;
    Kdc(k) = L.Kdc;
    Cn = pidtune(G{k}/Kdc(k), 'PI', 2*pi*fc, opt);
    Kpn(k) = Cn.Kp;
    Kin(k) = Cn.Ki;
end

grp = {1:6, 7:12};
gname = {'G2V','V2G'};
ns  = zeros(nm,1);
wg  = -inf(nm,1);
ref = zeros(nm,1);
for g = 1:2
    idx = grp{g};
    for r = idx
        st = false(numel(idx),1);
        gm = -inf(numel(idx),1);
        for j = 1:numel(idx)
            k = idx(j);
            R = loop_eval(G{k}, Kpn(r)/Kdc(k), Kin(r)/Kdc(k), P.Ts);
            st(j) = R.stable;
            gm(j) = R.gm;
        end
        ns(r) = sum(st);
        if any(st), wg(r) = min(gm(st)); end
    end
    [~, o] = sortrows([ns(idx) wg(idx)], [-1 -2]);
    ref(idx) = idx(o(1));
    fprintf('%s PI from %s:  Kp,n = %.4g   Ki,n = %.4g\n', gname{g}, M.label{idx(o(1))}, Kpn(idx(o(1))), Kin(idx(o(1))));
end
fprintf('\n');

S = cell(nm,1);
for k = 1:nm
    S{k} = loop_eval(G{k}, Kpn(ref(k))/Kdc(k), Kin(ref(k))/Kdc(k), P.Ts);
end
Kp  = Kpn(ref)./Kdc;
Ki  = Kin(ref)./Kdc;
fcx = cellfun(@(s) s.fc, S);
gmv = cellfun(@(s) s.gm, S);
pmv = cellfun(@(s) s.pm, S);
ncr = cellfun(@(s) s.nc, S);
rz  = cellfun(@(s) s.rz, S);
stb = cellfun(@(s) s.stable, S);

T0 = table(string(M.label(:)), Kpn, Kin, ns, wg, ...
    'VariableNames', {'Mode','Kp_n','Ki_n','StableModesWhenReference','WorstGainMargin_dB'});
T = table(string(M.label(:)), Kdc, Kp, Ki, fcx, gmv, pmv, ncr, rz, stb, ...
    'VariableNames', {'Mode','Kdc_A','Kp','Ki','Crossover_Hz','GainMargin_dB', ...
                      'PhaseMargin_deg','UnityCrossings','MaxClosedLoopPole','Stable'});

disp(T0);
disp(T);
writetable(T0, fullfile(rdir, 'pi_reference_scan.csv'));
writetable(T,  fullfile(rdir, 'pi_results.csv'));

function R = loop_eval(G, Kp, Ki, Ts)
C = pid(Kp, Ki, 'Ts', Ts, 'IFormula', 'ForwardEuler');
L = C*G;
[gm, pm, ~, wcp] = margin(L);
A = allmargin(L);
R.gm = 20*log10(gm);
R.pm = pm;
R.fc = wcp/(2*pi);
R.nc = numel(A.PMFrequency);
R.rz = max(abs(pole(feedback(L, 1))));
R.stable = R.rz < 1;
end