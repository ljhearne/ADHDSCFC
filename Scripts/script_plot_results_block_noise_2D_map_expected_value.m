% Load the file with the results from the model

clearvars
close all

%path = './';
path = '/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/';
load([path,'Results/Schaeffer214-Model/CTRL_block_evl_noise_25x25.mat'])

figure_handle = figure('Color','w','Position',[50 850 640 400]);

% Make axes for subplots
for kk=1:4
    ax(kk) = subplot(1, 4, kk);
    hold(ax(kk), 'on')
end


% Expected value of sigma in the hubs
esh = unique(SH);
% Expected value of sigma in the periphery
esp = unique(SP);

title_str = {'connections: all', 'hubs', 'feeders', 'local'};
midline_values  = {0.8914, 0.892, 0.8787, 0.9163}; 

rscfc_ctrl = [0.2708, 0.2394, 0.2318, 0.2933];
rscfc_adhd = [0.2571, 0.2350, 0.2100, 0.2838];
rscfc_percentage = ((rscfc_ctrl./  rscfc_adhd)-1);

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
    ax(kk).XTick = [0.4,0.6,0.8,1];
    ax(kk).YTick = [0.4,0.6,0.8,1];
    ax(kk).TickLength = [0.01,0.025];
    caxis(ax(kk), [0.65 0.95])
    axis(ax(kk), 'square')
    ax(kk).Title.String = title_str{kk};
    ax(kk).XLabel.String = 'E[\sigma_H]';
    ax(kk).YLabel.String = 'E[\sigma_P]';
    ch(kk) = colorbar(ax(kk));
    ch(kk).Location = 'SouthOutside';
    ch(kk).LineWidth = 1.5;    
end


% Plot identity line
line_color = {'k', 'k', 'k', 'k'};
for kk=1:4
    plot(ax(kk), esh, esp, 'color', line_color{kk}, 'linestyle', '-')
end
% Lines at 25% assymetry between hub and periphery
for kk=1:4
    %plot(ax(kk), esh, 1.25*esh, 'color', line_color{kk}, 'linestyle', '-')
    %plot(ax(kk), esh, 0.75*esh, 'color', line_color{kk}, 'linestyle', '-')

end
saveas(gcf,'/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/Results/K15/Schaefer214/Model_ScenarioI.svg');
%%  Plot the results for feeder edges only
% figure_handle_feeder = figure;
% ax_f = subplot(1,1,1);
% 
% ih(5) = imagesc(ax_f, esh, esp, r(3).map);
% set(ax_f, 'YDir', 'Normal')
% ax_f.XLim = [min(esh(:)) max(esh(:))];
% ax_f.YLim = [min(esp(:)) max(esp(:))];
% caxis(ax_f, [0.65 0.95])
% ax_f.Title.String = title_str{3};
% ax_f.XLabel.String = 'E[\sigma_H]';
% ax_f.YLabel.String = 'E[\sigma_P]';
% axis square
% hold(ax_f, 'on')
% plot(ax_f, esh, esp, 'k')
