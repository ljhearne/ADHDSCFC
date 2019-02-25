% Load one of these files and run the script
%load('./Results/Schaeffer214-Model/CTRL_block_var_noise_17x17_evals_adhd.mat')

%load('./Results/Schaeffer214-Model/CTRL_block_var_noise_16x16_evals_ctrl.mat')
load('./Results/Schaeffer214-Model/CTRL_block_var_noise_16x16_evals_adhd.mat')


figure_handle = figure;

% 
for kk=1:4
    ax(kk) = subplot(1, 4, kk);
end


% Get axis data and square to get variance info
% Variance in the hubs
vsh = unique(VSH).^2;
% Variance in the periphery
vsp = unique(VSP).^2;

title_str = {'all', 'h-h', 'feeder', 'p-p'};

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
    caxis(ax(kk), [min_r_val max_r_val])
    ax(kk).Title.String = title_str{kk};
    ax(kk).XLabel.String = 'Var[\sigma_H]';
    ax(kk).YLabel.String = 'Var[\sigma_P]';
end

% Colornar object
ch = colorbar;
ch.LineWidth = 1.5;
ch.AxisLocationMode = 'manual';
ch.Position = [0.9301 0.2057 0.0106 0.6971];

% Plot the results for feeder locations
figure_handle_feeder = figure;
ax_f = subplot(1,1,1);

ih(5) = imagesc(ax_f, vsh, vsp, imgaussfilt(r(3).map, 5));
set(ax_f, 'YDir', 'Normal')
ax_f.XLim = [min(vsh(:)) max(vsh(:))];
ax_f.YLim = [min(vsp(:)) max(vsp(:))];
caxis(ax_f, [min(r(3).map(:)) max(r(3).map(:))])
ax_f.Title.String = title_str{3};
ax_f.XLabel.String = 'Var[\sigma_H]';
ax_f.YLabel.String = 'Var[\sigma_P]';
axis square
hold(ax_f, 'on')
plot(ax_f, vsh, vsp, 'k')

%plot(ax_f, vhc, vpc, 'wx', 'markersize', 14)
plot(ax_f, vha, vpa, 'rx', 'markersize', 14)

rva = vha ./ vpa;
rvc = vhc ./ vpc;


%plot(ax_f, rva.*vsp, vsp, 'color', [1 0 0 0.1])
%plot(ax_f, rvc.*vsp, vsp, 'color', [1 1 1 0.1])


ch_f = colorbar;
ch_f.LineWidth = 1.5;
ch_f.AxisLocationMode = 'manual';
ch_f.Position = [0.9301 0.2057 0.0106 0.6971];
