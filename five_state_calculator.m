clear; close all; clc;
root = fileparts(mfilename('fullpath'));
fdir = fullfile(root, 'figures');
if ~isfolder(fdir), mkdir(fdir); end

[P, M] = dab_params();
pick = [1 7];
ns   = 40;
yl   = {'i_p (A)', 'i_{Li} (A)', 'i_{Lo} (A)', 'v_{Ci} (V)', 'v_{Co} (V)'};

f = figure('Color','w','Position',[60 40 1200 1000]);
tl = tiledlayout(f, 5, 2, 'TileSpacing','compact', 'Padding','compact', 'TileIndexing','columnmajor');

for c = 1:2
    k  = pick(c);
    D0 = [M.D(k,1); M.D(k,2); mod(M.D(k,3),2)];
    V0 = [P.V1; P.V2];
    [~, G, H] = dab_period_map(zeros(5,1), D0, V0, P);
    xs = (eye(5) - G)\H;
    [x1, ~, ~, ta, Xa] = dab_period_map(xs, D0, V0, P, ns);
    [~, ~, ~, tb, Xb] = dab_period_map(x1, D0, V0, P, ns);
    tt = [ta, ta(end) + tb(2:end)]*1e6;
    XX = [Xa, Xb(:,2:end)];
    fprintf('%s  D = [%.2f %.2f %.2f]\n', M.label{k}, M.D(k,:));
    for s = 1:5
        ax = nexttile(tl);
        plot(ax, tt, XX(s,:), '-', 'Color', [0 0 0.90], 'LineWidth', 1.4);
        xlim(ax, [0 2*P.Ts*1e6]);
        set(ax, 'FontName','Times New Roman', 'FontSize', 11, 'Box','on', 'TickDir','in', ...
                'LineWidth', 0.75, 'XColor','k', 'YColor','k');
        ax.Toolbar.Visible = 'off';
        if s == 1
            title(ax, sprintf('(%c)', 'a' + c - 1), 'FontWeight','normal', ...
                  'FontName','Times New Roman', 'FontSize', 12);
        end
        if s == 5, xlabel(ax, 't (\mus)'); end
        if c == 1, ylabel(ax, yl{s}); end
        fprintf('   %-11s mean %9.3f   min %9.3f   max %9.3f\n', yl{s}, mean(XX(s,:)), min(XX(s,:)), max(XX(s,:)));
    end
end

drawnow;
exportgraphics(f, fullfile(fdir, 'fig_states.png'), 'Resolution', 300);
fprintf('Saved %s\n', fullfile(fdir, 'fig_states.png'));