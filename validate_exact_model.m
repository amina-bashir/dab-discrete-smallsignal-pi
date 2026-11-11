clear; close all; clc;
root = fileparts(mfilename('fullpath'));
fdir = fullfile(root, 'figures');
if ~isfolder(fdir), mkdir(fdir); end

mdl = 'openloopDAB';
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

out = sim(mdl, 'StopTime', '0.01');

D0 = [D1; D2; mod(D3,2)];
V0 = [P.V1; P.V2];
[~, G, H] = dab_period_map(zeros(5,1), D0, V0, P);
xs = (eye(5) - G)\H;
[~, ~, ~, te, Xe] = dab_period_map(xs, D0, V0, P, 200);
te = te*1e6;

sig = {out.ip_sim, out.iLi_sim, out.iLo_sim, out.vCi_sim, out.vCo_sim};
yl  = {'i_p (A)', 'i_{Li} (A)', 'i_{Lo} (A)', 'v_{Ci} (V)', 'v_{Co} (V)'};

f = figure('Color','w','Position',[80 40 900 1100]);
tl = tiledlayout(f, 5, 1, 'TileSpacing','compact', 'Padding','compact');
fprintf('%s  D = [%.2f %.2f %.2f]\n', M.label{k}, D1, D2, D3);
for q = 1:5
    t = sig{q}.Time(:);
    y = squeeze(sig{q}.Data);
    y = y(:);
    np = floor(t(end)/P.Ts + 1e-9);
    t0 = (np - 1)*P.Ts;
    m  = t >= t0 - 1e-12 & t <= t0 + P.Ts + 1e-12;
    ts = (t(m) - t0)*1e6;
    ys = y(m);
    [tu, iu] = unique(ts);
    yi  = interp1(tu, ys(iu), te, 'linear', 'extrap');
    ye  = Xe(q,:);
    err = 100*sqrt(mean((yi - ye).^2))/max(abs(ye));
    fprintf('   %-10s  mean sim %9.3f   mean model %9.3f   rms difference %.3f %% of peak\n', ...
            yl{q}, mean(yi), mean(ye), err);

    ax = nexttile(tl);
    h1 = plot(ax, ts, ys, '-', 'Color', [0.85 0.60 0], 'LineWidth', 1.8);
    hold(ax, 'on');
    h2 = plot(ax, te, ye, '--', 'Color', [0 0 0.90], 'LineWidth', 1.3);
    xlim(ax, [0 P.Ts*1e6]);
    set(ax, 'FontName','Times New Roman', 'FontSize', 11, 'Box','on', 'TickDir','in', ...
            'LineWidth', 0.75, 'XColor','k', 'YColor','k');
    ax.Toolbar.Visible = 'off';
    title(ax, sprintf('(%c)', 'a' + q - 1), 'FontWeight','normal', ...
          'FontName','Times New Roman', 'FontSize', 12);
    if q == 5, xlabel(ax, 't (\mus)'); end
    ylabel(ax, yl{q});
end
lg = legend([h1 h2], {'Simulink', 'Exact discrete-time model'}, 'Orientation','horizontal', ...
            'FontName','Times New Roman', 'FontSize', 11);
lg.Layout.Tile = 'north';

drawnow;
exportgraphics(f, fullfile(fdir, 'fig_simulink_validation.png'), 'Resolution', 300);
fprintf('Saved %s\n', fullfile(fdir, 'fig_simulink_validation.png'));