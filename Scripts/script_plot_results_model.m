% This script plot the figures using the correlation values calculated 
% from eSC-aFC when all the nodes have noise variability.  
close all
clearvars

DocsPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/';
pathtofiles = [DocsPath,'Results' filesep 'Schaeffer214-Model' filesep];
resultsdir = [DocsPath,'Results/K15/Schaefer214' filesep];
filenames = {'CTRL_variable_noise_all_nodes_seed_all', ...
             'CTRL_variable_noise_hub_nodes_seed_all', ...
             'CTRL_variable_noise_periphery_nodes_seed_all'};

for ii=1:3
    load([pathtofiles filenames{ii}])
end
% median over seeds, and then median over subjects
ctrl_all_std_values       = median(median(all_std_values, 3), 1);         
ctrl_periphery_std_values = median(median(periphery_std_values, 3), 1);
ctrl_hubs_std_values      = median(median(hubs_std_values, 3), 1);

%% Plot in format for illustrator
figure('Color','w','Position',[450 450 400 175]);

subplot(1,3,1)
% HUB edges
plot(ctrl_all_std_values,median(median(r_ctrl_esc_afc.hub_edges, 3), 1),...
    'color', [0.5 0.5 0.5], 'linestyle', '-','linewidth',2); hold on
plot(ctrl_hubs_std_values,median(median(r_ctrl_esc_afc_hubs.hub_edges, 3), 1),...
    'color', [1 0 0], 'linestyle', '-','linewidth',2); hold on
plot(ctrl_periphery_std_values,median(median(r_ctrl_esc_afc_periphery.hub_edges, 3), 1),...
    'color', [1 1 0], 'linestyle', '-','linewidth',2); hold on

set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
set(gca,'YLim',[0.2 1]);
set(gca,'XLim',[0 0.5]);
set(gca,'XTick',0.1:0.1:0.5);
%xlabel('Std Dev [\sigma_{i}]'); % This is \bar{\nu} in the text
%ylabel(y_label);

subplot(1,3,2)
% FEEDER edges
plot(ctrl_all_std_values,median(median(r_ctrl_esc_afc.feed_edges, 3), 1), ...
    'color', [0.5 0.5 0.5], 'linestyle', '-','linewidth',2); hold on
plot(ctrl_hubs_std_values,median(median(r_ctrl_esc_afc_hubs.feed_edges, 3), 1), ...
    'color', [1 0 0], 'linestyle', '-','linewidth',2); hold on
plot(ctrl_periphery_std_values,median(median(r_ctrl_esc_afc_periphery.feed_edges, 3), 1), ...
    'color', [1 1 0], 'linestyle', '-','linewidth',2); hold on

subplot(1,3,3)
% LOCAL edges
plot(ctrl_all_std_values,median(median(r_ctrl_esc_afc.periphery_edges, 3), 1), ...
    'color', [0.5 0.5 0.5], 'linestyle', '-','linewidth',2); hold on
plot(ctrl_hubs_std_values,median(median(r_ctrl_esc_afc_hubs.periphery_edges, 3), 1), ...
    'color', [1 0 0], 'linestyle', '-','linewidth',2); hold on
plot(ctrl_periphery_std_values,median(median(r_ctrl_esc_afc_periphery.periphery_edges, 3), 1), ...
    'color', [1 1 0], 'linestyle', '-','linewidth',2); hold on

set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
set(gca,'YLim',[0.2 1]);
set(gca,'XLim',[0 0.5]);
set(gca,'XTick',0.1:0.1:0.5);
%xlabel('Std Dev [\sigma_{i}]'); % This is \bar{\nu} in the text
%ylabel(y_label);

saveas(gcf,[resultsdir,'Model_results.svg']);

%% plot the variability across simulations...
figure('Color','w','Position',[850 450 400 175]);

subplot(1,3,1)
% HUB edges
model_plot(ctrl_all_std_values,r_ctrl_esc_afc.hub_edges,[.5,.5,.5])
model_plot(ctrl_hubs_std_values,r_ctrl_esc_afc_hubs.hub_edges,[1,0,0])
model_plot(ctrl_periphery_std_values,r_ctrl_esc_afc_periphery.hub_edges,[1,.78,0])

set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
set(gca,'YLim',[0.2 1]);
set(gca,'XLim',[0 0.5]);
set(gca,'XTick',0.1:0.2:0.5);
%xlabel('Std Dev [\sigma_{i}]'); % This is \bar{\nu} in the text
%ylabel(y_label);

subplot(1,3,2)
% FEEDER edges
model_plot(ctrl_all_std_values,r_ctrl_esc_afc.feed_edges,[.5,.5,.5])
model_plot(ctrl_hubs_std_values,r_ctrl_esc_afc_hubs.feed_edges,[1,0,0])
model_plot(ctrl_periphery_std_values,r_ctrl_esc_afc_periphery.feed_edges,[1,.78,0])
set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
set(gca,'YLim',[0.2 1]);
set(gca,'XLim',[0 0.5]);
set(gca,'XTick',0.1:0.2:0.5);

subplot(1,3,3)
% LOCAL edges
model_plot(ctrl_all_std_values,r_ctrl_esc_afc.periphery_edges,[.5,.5,.5])
model_plot(ctrl_hubs_std_values,r_ctrl_esc_afc_hubs.periphery_edges,[1,0,0])
model_plot(ctrl_periphery_std_values,r_ctrl_esc_afc_periphery.periphery_edges,[1,.78,0])
set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
set(gca,'YLim',[0.2 1]);
set(gca,'XLim',[0 0.5]);
set(gca,'XTick',0.1:0.2:0.5);


%xlabel('Std Dev [\sigma_{i}]'); % This is \bar{\nu} in the text
%ylabel(y_label);
saveas(gcf,[resultsdir,'Model_resultsPermutations.svg']);

function model_plot(x,y,color)
for i = 1:size(y,3)
    h = plot(x,median(y(:,:,i),1),...
    'color', color, 'linestyle', '-','linewidth',1); hold on
    h.Color(4) = 0.1;  % 50% transparent
end
plot(x,median(median(y, 3), 1),...
    'color', color, 'linestyle', '-','linewidth',2); hold on
end