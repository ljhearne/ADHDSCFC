% Load modelling results
load('./Results/Schaeffer214-Model/CTRL_block_var_noise_16x16_evals_ctrl.mat')


figure_handle = figure;

% 
for kk=1:4
    ax(kk) = subplot(1, 4, kk);
    hold(ax(kk), 'on')
end


% Get axis data and square to get variance info
% Variance in the hubs
vsh = unique(VSH).^2;
% Variance in the periphery
vsp = unique(VSP).^2;

title_str = {'connections: all', 'connections: h-h', 'connections: feeder', 'connections: p-p'};

% Get results for all type of connections
% All connections
r(1).map = reshape(median(r_ctrl_esc_afc.all_edges), size(VSH));
% Hub-Hub connections
r(2).map = reshape(median(r_ctrl_esc_afc.hub_edges), size(VSH));
% Feeder
r(3).map = reshape(median(r_ctrl_esc_afc.feed_edges), size(VSH));
% Periphery
r(4).map = reshape(median(r_ctrl_esc_afc.periphery_edges), size(VSH));

min_r_val = [];
max_r_val = [];

for kk=1:4
    min_r_val = [min_r_val min(r(kk).map)];
    max_r_val = [max_r_val max(r(kk).map)];
end

min_r_val = min(min_r_val);
max_r_val = max(max_r_val);


for kk=1:4
    ih(kk) = imagesc(ax(kk), vsh, vsp, r(kk).map);
    set(ax(kk), 'YDir', 'Normal')
    ax(kk).XLim = [min(vsh(:)) max(vsh(:))];
    ax(kk).YLim = [min(vsp(:)) max(vsp(:))];
    caxis(ax(kk), [0.2 0.95])
    axis(ax(kk), 'square')
    ax(kk).Title.String = title_str{kk};
    ax(kk).XLabel.String = 'Var[\sigma_H]';
    ax(kk).YLabel.String = 'Var[\sigma_P]';
    ch(kk) = colorbar(ax(kk));
    ch(kk).Location = 'SouthOutside';
    ch(kk).LineWidth = 1.5;    
end

% Plot identity line
line_color = {'k', 'w', 'k', 'w'};
for kk=1:4
    plot(ax(kk), vsh, vsp, 'color', line_color{kk}, 'linestyle', '-')
end

%% Plot the results for feeder locations
figure_handle_feeder = figure;
ax_f = subplot(1,1,1);

ih(5) = imagesc(ax_f, vsh, vsp, r(3).map);
set(ax_f, 'YDir', 'Normal')
ax_f.XLim = [min(vsh(:)) max(vsh(:))];
ax_f.YLim = [min(vsp(:)) max(vsp(:))];
caxis(ax_f, [0.2 0.95])
ax_f.Title.String = title_str{3};
ax_f.XLabel.String = 'Var[\sigma_H]';
ax_f.YLabel.String = 'Var[\sigma_P]';
axis square
hold(ax_f, 'on')
plot(ax_f, vsh, vsp, 'k')
ch_f = colorbar;
