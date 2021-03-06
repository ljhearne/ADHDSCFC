% Plot 2D var maps for all connections
clearvars
close all

%path = './';
path = '/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/';
load([path,'Results/Schaeffer214-Model/CTRL_block_var_noise_16x16_hub_periphery_10_percent.mat'])

figure_handle = figure('Color','w','Position',[50 850 640 400]);
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

cmap = parula(65);
for kk=1:4
    dummy = imgaussfilt(r(kk).map, 1);
    ih(kk) = imagesc(ax(kk), vsh, vsp, dummy);

    set(ax(kk), 'YDir', 'Normal')
    ax(kk).XLim = [min(vsh(:)) max(vsh(:))];
    ax(kk).YLim = [min(vsp(:)) max(vsp(:))];
    ax(kk).XTick = [0.05,0.15,0.25];
    ax(kk).YTick = [0.05,0.15,0.25];
    caxis(ax(kk), [0.15 1])
    axis(ax(kk), 'square')
    ax(kk).Title.String = title_str{kk};
    ax(kk).XLabel.String = 'Var[\sigma_H]';
    ax(kk).YLabel.String = 'Var[\sigma_P]';
    ch(kk) = colorbar(ax(kk));
    ch(kk).Location = 'SouthOutside';
    ch(kk).LineWidth = 1.5;
    ax(kk).Colormap = cmap;
end

% Plot identity line
line_color = {'k', 'k', 'k', 'k'};
for kk=1:4
    plot(ax(kk), vsh, vsp, 'color', line_color{kk}, 'linestyle', '-')
end
saveas(gcf,'/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/Results/K15/Schaefer214/Model_ScenarioII.svg');
% %% Plot the results for feeder connections
% figure_handle_feeder = figure;
% ax_f = subplot(1,1,1);
% 
% ih(5) = imagesc(ax_f, vsh, vsp, imgaussfilt(r(3).map, 1));
% set(ax_f, 'YDir', 'Normal')
% ax_f.XLim = [min(vsh(:)) max(vsh(:))];
% ax_f.YLim = [min(vsp(:)) max(vsp(:))];
% caxis(ax_f, [0.15 1])
% ax_f.Title.String = title_str{3};
% ax_f.XLabel.String = 'Var[\sigma_H]';
% ax_f.YLabel.String = 'Var[\sigma_P]';
% axis square
% hold(ax_f, 'on')
% plot(ax_f, vsh, vsp, 'k')
% ch_f = colorbar;