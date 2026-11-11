clear; close all; clc;
root = fileparts(mfilename('fullpath'));
fdir = fullfile(root, 'figures');
if ~isfolder(fdir), mkdir(fdir); end

[P, M] = dab_params();
nm = numel(M.name);
ns = 40;

f = figure('Color','w','Position',[60 60 1500 850]);
tl = tiledlayout(f, 3, 4, 'TileSpacing','compact', 'Padding','compact');

for k = 1:nm
    D0 = [M.D(k,1); M.D(k,2); mod(M.D(k,3),2)];
    V0 = [P.V1; P.V2];
    [~, G, H] = dab_period_map(zeros(5,1), D0, V0, P);
    xs = (eye(5) - G)\H;
    [x1, ~, ~, ta, Xa] = dab_period_map(xs, D0, V0, P, ns);
    [~, ~, ~, tb, Xb] = dab_period_map(x1, D0, V0, P, ns);
    tt = [ta, ta(end) + tb(2:end)]*1e6;
    ii = [Xa(1,:), Xb(1,2:end)];

    ax = nexttile(tl);
    plot(ax, tt, ii, '-', 'Color', [0 0 0.90], 'LineWidth', 1.4);
    xlim(ax, [0 2*P.Ts*1e6]);
    fmt(ax, k, 't (\mus)', 'i_p (A)');
    fprintf('%-14s  peak %8.2f A   rms %8.2f A\n', M.label{k}, max(abs(ii)), sqrt(mean(ii.^2)));
end

save_png(f, fullfile(fdir, 'fig_link_current.png'));

function fmt(ax, k, xl, yl)
set(ax, 'FontName','Times New Roman', 'FontSize', 11, 'Box','on', 'TickDir','in', ...
        'LineWidth', 0.75, 'XColor','k', 'YColor','k');
title(ax, sprintf('(%c)', 'a' + k - 1), 'FontWeight','normal', ...
      'FontName','Times New Roman', 'FontSize', 12);
if k > 8, xlabel(ax, xl); end
if mod(k-1, 4) == 0, ylabel(ax, yl); end
end

function save_png(f, fn)
try
    exportgraphics(f, fn, 'Resolution', 300);
catch
    print(f, fn, '-dpng', '-r300');
end
end
