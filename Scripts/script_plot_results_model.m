% This script plot the figures using the correlation values calculated 
% from eSC-aFC when all the nodes have noise variability.  

% Assuming you're at the top level of the ADHDSCFC folder
pathtofiles = ['Results' filesep 'Schaeffer214-Model' filesep];
filenames = {'CTRL_variable_noise_all_nodes_seed_all', ...
             'CTRL_variable_noise_hub_nodes_seed_all', ...
             'CTRL_variable_noise_periphery_nodes_seed_all'};

for ii=1:3
    load([pathtofiles filenames{ii}])
end
%%
% median over seeds, and then median over subjects
ctrl_all_std_values       = median(median(all_std_values, 3), 1);         
ctrl_periphery_std_values = median(median(periphery_std_values, 3), 1);
ctrl_hubs_std_values      = median(median(hubs_std_values, 3), 1);

%% Plot Fig. 6 - Supplementary material
figure_handle_a = figure;

% Each subplot will plot the correlation curves based on
ax(1) = subplot(2, 2, 1); % all edges
ax(2) = subplot(2, 2, 2); % hub edges 
ax(3) = subplot(2, 2, 3); % feeder edges
ax(4) = subplot(2, 2, 4); % periphery edges


titles_subplots = {'r - all connections', ...
                   'r - hub connections', ...
                   'r - feeder connections', ...
                   'r - periphery connections'};

for ii=1:length(ax)
   hold(ax(ii), 'on')
   ax(ii).Title.String = titles_subplots{ii};
end

y_label = {'r_{eSC-aFC}'};

% ALL edges
h(1, 1) = plot(ax(1), ctrl_all_std_values,       median(median(r_ctrl_esc_afc.all_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-');
h(1, 2) = plot(ax(1), ctrl_hubs_std_values,      median(median(r_ctrl_esc_afc_hubs.all_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '--');
h(1, 3) = plot(ax(1), ctrl_periphery_std_values, median(median(r_ctrl_esc_afc_periphery.all_edges, 3),1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-.');
% HUB edges
h(2, 1) = plot(ax(2), ctrl_all_std_values,       median(median(r_ctrl_esc_afc.hub_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-');
h(2, 2) = plot(ax(2), ctrl_hubs_std_values,      median(median(r_ctrl_esc_afc_hubs.hub_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '--');
h(2, 3) = plot(ax(2), ctrl_periphery_std_values, median(median(r_ctrl_esc_afc_periphery.hub_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-.');
% FEEDER edges
h(3, 1) = plot(ax(3), ctrl_all_std_values,       median(median(r_ctrl_esc_afc.feed_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-');
h(3, 2) = plot(ax(3), ctrl_hubs_std_values,      median(median(r_ctrl_esc_afc_hubs.feed_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '--');
h(3, 3) = plot(ax(3), ctrl_periphery_std_values, median(median(r_ctrl_esc_afc_periphery.feed_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-.');

% PERIPHERY edges
h(4, 1) = plot(ax(4), ctrl_all_std_values,       median(median(r_ctrl_esc_afc.periphery_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-');
h(4, 2) = plot(ax(4), ctrl_hubs_std_values,      median(median(r_ctrl_esc_afc_hubs.periphery_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '--');
h(4, 3) = plot(ax(4), ctrl_periphery_std_values, median(median(r_ctrl_esc_afc_periphery.periphery_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-.');


for ii=1:length(ax)
    ax(ii).YLim = [0.45 1];
    ax(ii).XLim = [0 0.5];
    ax(ii).XLabel.String = 'Std Dev [\sigma_{i}]'; % This is \bar{\nu} in the text
    ax(ii).YLabel.String = y_label;
end


%% Plot CTRL only for hub and feeder edges 
figure_handle_b = figure;

% Each subplot will plot the correlation curves based on
bx(1) = subplot(1, 2, 1); % hub    edges
bx(2) = subplot(1, 2, 2); % feeder edges 


titles_subplots = {'r - hub connections', ...
                   'r - feeder connections'};

for ii=1:length(bx)
   hold(bx(ii), 'on')
   bx(ii).Title.String = titles_subplots{ii};
end


% HUB edges
hb(1, 1) = plot(bx(1), ctrl_all_std_values,       median(median(r_ctrl_esc_afc.hub_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-');
hb(1, 2) = plot(bx(1), ctrl_hubs_std_values,      median(median(r_ctrl_esc_afc_hubs.hub_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '--');
hb(1, 3) = plot(bx(1), ctrl_periphery_std_values, median(median(r_ctrl_esc_afc_periphery.hub_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-.');
% FEEDER edges
hb(1, 1) = plot(bx(2), ctrl_all_std_values,       median(median(r_ctrl_esc_afc.feed_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-');
hb(1, 2) = plot(bx(2), ctrl_hubs_std_values,      median(median(r_ctrl_esc_afc_hubs.feed_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '--');
hb(1, 3) = plot(bx(2), ctrl_periphery_std_values, median(median(r_ctrl_esc_afc_periphery.feed_edges, 3), 1), ...
                                                 'color', [0.5 0.5 0.5], 'linestyle', '-.');

y_label = {'r_{eSC-aFC}'};
for ii=1:length(bx)
    bx(ii).YLim = [0.2 1];
    bx(ii).XLim = [0 0.5];
    bx(ii).XLabel.String = 'Std Dev [\sigma_{i}]'; % This is \bar{\nu} in the text
    bx(ii).YLabel.String = y_label;
end

