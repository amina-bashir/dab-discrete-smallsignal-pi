clear; close all; clc;
root = fileparts(mfilename('fullpath'));
fdir = fullfile(root, 'figures');
rdir = fullfile(root, 'results');
if ~isfolder(fdir), mkdir(fdir); end
if ~isfolder(rdir), mkdir(rdir); end

[P, M] = dab_params();
nm = numel(M.name);
N0 = 50;
Nc = 400;
dD = 1e-3;
fq = logspace(1, log10(P.fs/2), 600);
cX = [0.85 0.60 0];
cL = [0 0 0.90];
cG = [0 0.30 0.85];
cV = [0.85 0.10 0.10];

f1 = figure('Color','w','Position',[40 40 1500 900]);
t1 = tiledlayout(f1, 3, 4, 'TileSpacing','compact', 'Padding','compact');
f2 = figure('Color','w','Position',[60 60 1500 900]);
t2 = tiledlayout(f2, 3, 4, 'TileSpacing','compact', 'Padding','compact');

Kdc = zeros(nm,1);
rho = zeros(nm,1);
err = zeros(nm,1);
pz  = cell(nm,2);

for k = 1:nm
    L = dab_linearize(M.D(k,:), P);
    Kdc(k) = L.Kdc;
    rho(k) = max(abs(eig(L.A)));
    pz{k,1} = pole(L.G);
    pz{k,2} = tzero(L.G);

    xn = L.xs;
    xl = zeros(5,1);
    yn = zeros(1,Nc);
    yl = zeros(1,Nc);
    for j = 1:Nc
        u = dD*(j >= N0);
        xl = L.A*xl + L.B(:,3)*u;
        yl(j) = L.xs(3) + xl(3);
        xn = dab_period_map(xn, L.D0 + [0;0;u], L.V0, P);
        yn(j) = xn(3);
    end
    err(k) = 100*norm(yn(N0:end) - yl(N0:end))/norm(yn(N0:end) - yn(N0-1));

    tm = (0:Nc-1)*P.Ts*1e3;
    ax = nexttile(t1);
    h1 = plot(ax, tm, yn, '-', 'Color', cX, 'LineWidth', 1.8);
    hold(ax, 'on');
    h2 = plot(ax, tm, yl, '--', 'Color', cL, 'LineWidth', 1.3);
    xlim(ax, [N0-15, N0+150]*P.Ts*1e3);
    fmt(ax, k, 't (ms)', 'i_{Lo} (A)');

    mag = squeeze(bode(L.G, 2*pi*fq));
    ax = nexttile(t2);
    semilogx(ax, fq, 20*log10(mag), '-', 'Color', cL, 'LineWidth', 1.4);
    xlim(ax, [fq(1) fq(end)]);
    fmt(ax, k, 'f (Hz)', '|G_{i_{Lo},D_3}| (dB)');
end

lg = legend([h1 h2], {'Exact discrete-time model', 'Small-signal model'}, ...
            'Orientation','horizontal', 'FontName','Times New Roman', 'FontSize', 11);
lg.Layout.Tile = 'north';
save_png(f1, fullfile(fdir, 'fig_step_validation.png'));
save_png(f2, fullfile(fdir, 'fig_bode.png'));

f3 = figure('Color','w','Position',[80 80 1300 520]);
t3 = tiledlayout(f3, 1, 2, 'TileSpacing','compact', 'Padding','compact');
th = linspace(0, 2*pi, 400);
mk = {'+','x'; 'o','o'};
allz = vertcat(pz{:,2});
for q = 1:2
    ax = nexttile(t3);
    hold(ax, 'on');
    plot(ax, cos(th), sin(th), 'k:', 'LineWidth', 0.8);
    hh = gobjects(nm,1);
    for k = 1:nm
        v = pz{k,q};
        if k <= 6, c = cG; d = 1; else, c = cV; d = 2; end
        hh(k) = plot(ax, real(v), imag(v), mk{q,d}, 'Color', c, 'MarkerSize', 9, 'LineWidth', 1.3);
    end
    daspect(ax, [1 1 1]);
    if q == 1
        xlim(ax, [-1.2 1.2]);
        ylim(ax, [-1.2 1.2]);
    else
        xlim(ax, [floor(min(real(allz))) - 0.5, 1.5]);
        ylim(ax, [-1.5 1.5]);
    end
    fmt(ax, q, 'Re(z)', 'Im(z)');
    xlabel(ax, 'Re(z)');
    ylabel(ax, 'Im(z)');
    if q == 1
        legend(ax, [hh(1) hh(7)], {'G2V','V2G'}, 'Location','southwest', ...
               'FontName','Times New Roman', 'FontSize', 11);
    end
end
save_png(f3, fullfile(fdir, 'fig_poles_zeros.png'));

T = table(string(M.label(:)), M.D(:,1), M.D(:,2), M.D(:,3), Kdc, rho, err, ...
    'VariableNames', {'Mode','D1','D2','D3','Kdc_A','MaxFloquetMultiplier','StepErrorPct'});
disp(T);
writetable(T, fullfile(rdir, 'smallsignal_summary.csv'));

function fmt(ax, k, xl, yl)
set(ax, 'FontName','Times New Roman', 'FontSize', 11, 'Box','on', 'TickDir','in', ...
        'LineWidth', 0.75, 'XColor','k', 'YColor','k');
ax.Toolbar.Visible = 'off';
title(ax, sprintf('(%c)', 'a' + k - 1), 'FontWeight','normal', ...
      'FontName','Times New Roman', 'FontSize', 12);
if k > 8, xlabel(ax, xl); end
if mod(k-1, 4) == 0, ylabel(ax, yl); end
end

function save_png(f, fn)
if ~isgraphics(f)
    fprintf('Not saved, figure window was closed: %s\n', fn);
    return
end
drawnow;
exportgraphics(f, fn, 'Resolution', 300);
fprintf('Saved %s\n', fn);
end