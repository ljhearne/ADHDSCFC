% Load 2D maps Var[\sigma_P] vs Var[\sigma_H]
clearvars
close all

addpath('Functions/');

%path = './';
path = '/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/';

load([path,'Results/Schaeffer214-Model/CTRL_block_var_noise_16x16_hub_periphery_00_percent.mat'])
map = reshape(median(r_ctrl_esc_afc.feed_edges), size(VSH));
feeder_map(1).map = map';
load([path,'Results/Schaeffer214-Model/CTRL_block_var_noise_16x16_hub_periphery_10_percent.mat'])
map = reshape(median(r_ctrl_esc_afc.feed_edges), size(VSH));
feeder_map(2).map = map';
load([path,'Results/Schaeffer214-Model/CTRL_block_var_noise_16x16_hub_periphery_20_percent.mat'])
map = reshape(median(r_ctrl_esc_afc.feed_edges), size(VSH));
feeder_map(3).map = map';
load([path,'Results/Schaeffer214-Model/CTRL_block_var_noise_16x16_hub_periphery_50_percent.mat'])
map = reshape(median(r_ctrl_esc_afc.feed_edges), size(VSH));
feeder_map(4).map = map';

figure_handle = figure('Color','w','Position',[50 850 800 500]);

% 
for kk=1:8
    ax(kk) = subplot(2, 4, kk);
    hold(ax(kk), 'on')
end


% Get axis data and square to get variance info
% Variance in the hubs
vsh = unique(VSH).^2;
% Variance in the periphery
vsp = unique(VSP).^2;

title_str = {'E[\sigma_H] = E[\sigma_P]', 'E[\sigma_H] < E[\sigma_P] (10%)', 'E[\sigma_H] < E[\sigma_P] (20%)', 'E[\sigma_H] < E[\sigma_P] (50%)'};
midline_values  = {0.8914, 0.892, 0.8787, 0.9163}; 

% Mean values from empirical data
rscfc_ctrl = [0.2708, 0.2394, 0.2318, 0.2933];
rscfc_adhd = [0.2571, 0.2350, 0.2100, 0.2838];
rscfc_percentage = ((rscfc_ctrl./  rscfc_adhd)-1);

cmap_corr = parula(64);
for kk=1:4
    dummy = imgaussfilt(feeder_map(kk).map, 1);
    ih(kk) = imagesc(ax(kk), vsh, vsp, dummy);
    
    set(ax(kk), 'YDir', 'Normal')
    ax(kk).XLim = [min(vsh(:)) max(vsh(:))];
    ax(kk).YLim = [min(vsp(:)) max(vsp(:))];
    caxis([0.15 1])
    axis(ax(kk), 'square')
    ax(kk).Title.String = title_str{kk};
    ax(kk).XLabel.String = 'Var[\sigma_H]';
    ax(kk).YLabel.String = 'Var[\sigma_P]';
    ch(kk) = colorbar(ax(kk));
    ch(kk).Location = 'SouthOutside';
    ch(kk).LineWidth = 1.5;
    ax(kk).Colormap = cmap_corr;
    ax(kk).FontSize=12;
    ax(kk).FontName='Helvetica';
end

% Plot identity line
line_color = {'k', 'k', 'k', 'k'};
for kk=1:4
    plot(ax(kk), vsh, vsp, 'color', line_color{kk}, 'linestyle', '-')
end
delete(ax(5))
%% Plot the difference maps

cmap_diff = bluered(65);
for kk=6:8
    dummy = imgaussfilt((feeder_map(kk-4).map./feeder_map(1).map)-1, 1);
    ih(kk) = imagesc(ax(kk), vsh, vsp, dummy*100);
    set(ax(kk), 'YDir', 'Normal')
    ax(kk).XLim = [min(vsh(:)) max(vsh(:))];
    ax(kk).YLim = [min(vsp(:)) max(vsp(:))];
    axis(ax(kk), 'square')        
    caxis(ax(kk), [-50 50])
    ax(kk).XLabel.String = 'Var[\sigma_H]';
    ax(kk).YLabel.String = 'Var[\sigma_P]';
    ch(kk) = colorbar(ax(kk));
    ch(kk).Location = 'SouthOutside';
    ch(kk).LineWidth = 1.5;
    ax(kk).Colormap = cmap_diff;
    plot(ax(kk), vsh, vsp, 'color', line_color{kk-4}, 'linestyle', '-')
end


