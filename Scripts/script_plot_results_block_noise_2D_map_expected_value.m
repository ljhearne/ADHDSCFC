% Load the file
load('./Results/Schaeffer214-Model/CTRL_block_evl_noise_25x25.mat')

figure_handle = figure;

% 
for kk=1:4
    ax(kk) = subplot(1, 4, kk);
    hold(ax(kk), 'on')
end


% Expected value of sigma in the hubs
esh = unique(SH);
% Expected value of sigma in the periphery
esp = unique(SP);

title_str = {'all', 'h-h', 'feeder', 'p-p'};

% Get results for all type of connections
% All connections
r(1).map = reshape(median(r_ctrl_esc_afc.all_edges), size(SH));
% Hub-Hub connections
r(2).map = reshape(median(r_ctrl_esc_afc.hub_edges), size(SH));
% Feeder
r(3).map = reshape(median(r_ctrl_esc_afc.feed_edges), size(SH));
% Periphery
r(4).map = reshape(median(r_ctrl_esc_afc.periphery_edges), size(SH));

min_r_val = [];
max_r_val = [];

for kk=1:4
    min_r_val = [min_r_val min(r(kk).map)];
    max_r_val = [max_r_val max(r(kk).map)];
end

min_r_val = min(min_r_val);
max_r_val = max(max_r_val);


for kk=1:4
    ih(kk) = imagesc(ax(kk), esh, esp, r(kk).map);
    set(ax(kk), 'YDir', 'Normal')
    ax(kk).XLim = [min(esh(:)) max(esh(:))];
    ax(kk).YLim = [min(esp(:)) max(esp(:))];
    caxis(ax(kk), [min(r(kk).map(:)) max(r(kk).map(:))])
    axis(ax(kk), 'square')
    ax(kk).Title.String = title_str{kk};
    ax(kk).XLabel.String = 'E[\sigma_H]';
    ax(kk).YLabel.String = 'E[\sigma_P]';
    ch(kk) = colorbar(ax(kk));
    ch(kk).Location = 'SouthOutside';
    ch(kk).LineWidth = 1.5;    
end


% Plot addittional things

line_color = {'k', 'w', 'k', 'w'};
for kk=1:4
    plot(ax(kk), esh, esp, 'color', line_color{kk}, 'linestyle', '-')
    plot(ax(kk), esh, 1.02*esh, 'color', line_color{kk}, 'linestyle', '--') %CTRL
    plot(ax(kk), esh, 1.1*esh, 'color', line_color{kk}, 'linestyle', ':')   % ADHD
end


% assume value for ctrl var(sigma_h) and var(sigma_p)
ev_sh_c = 0.5;
ev_sp_c = ev_sh_c*1.01;
ev_sh_a = 0.97*ev_sh_c;
ev_sp_a = 1.08*ev_sp_c;

% assume value for ctrl var(sigma_h) and var(sigma_p)
var_sh_c = 0.04;
var_sp_c = var_sh_c * 0.75;
var_sh_a = var_sh_c + var_sh_c * 0.03;
var_sp_a = var_sp_c + var_sp_c * 0.3;

[x_c, y_c] = calculate_ellipse_line(ev_sh_c, ev_sp_c, var_sh_c, var_sp_c, 0);
[x_a, y_a] = calculate_ellipse_line(ev_sh_a, ev_sp_a, var_sh_a, var_sp_a, 0);


plot(ax(3), ev_sh_c, ev_sp_c, 'gx', 'markersize', 14)
plot(ax(3), ev_sh_a, ev_sp_a, 'rx', 'markersize', 14)
plot(ax(3), x_c, y_c, 'g')
plot(ax(3), x_a, y_a, 'r')



ch_f = colorbar;
ch_f.LineWidth = 1.5;
ch_f.AxisLocationMode = 'manual';
ch_f.Position = [0.9301 0.2057 0.0106 0.6971];
