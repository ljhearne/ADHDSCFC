% Load the file with the results from the model
load('./Results/Schaeffer214-Model/CTRL_block_evl_noise_25x25.mat')

figure_handle = figure;

% Make axes for subplots
for kk=1:4
    ax(kk) = subplot(1, 4, kk);
    hold(ax(kk), 'on')
end


% Expected value of sigma in the hubs
esh = unique(SH);
% Expected value of sigma in the periphery
esp = unique(SP);

title_str = {'connections: all', 'connections: h-h', 'connections: feeder', 'connections: p-p'};

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
line_color = {'k', 'w', 'k', 'w'};
for kk=1:4
    plot(ax(kk), esh, esp, 'color', line_color{kk}, 'linestyle', '-')
end

%% 
line_color = {'k', 'w', 'k', 'w'};
for kk=1:4
    plot(ax(kk), esh, 1.15*esh, 'color', line_color{kk}, 'linestyle', '-')
    plot(ax(kk), esh, 0.85*esh, 'color', line_color{kk}, 'linestyle', '-')

end
%%  Plot the results for feeder edges only
figure_handle_feeder = figure;
ax_f = subplot(1,1,1);

ih(5) = imagesc(ax_f, esh, esp, r(3).map);
set(ax_f, 'YDir', 'Normal')
ax_f.XLim = [min(esh(:)) max(esh(:))];
ax_f.YLim = [min(esp(:)) max(esp(:))];
caxis(ax_f, [0.65 0.95])
ax_f.Title.String = title_str{3};
ax_f.XLabel.String = 'E[\sigma_H]';
ax_f.YLabel.String = 'E[\sigma_P]';
axis square
hold(ax_f, 'on')
plot(ax_f, esh, esp, 'k')

%% Plot the assymetry mapsplot(ax_h, esh, 1.15*esh, 'color', 'k', 'linestyle', '--')

figure_handle_hub_asym = figure;
ax_h = subplot(1,1,1);
hold(ax_h, 'on')
hub_map = r(4).map';

% Calculate differences
tri_low = tril(hub_map);
tri_up = triu(hub_map)';
res = tri_up - tri_low;

im_handle = imagesc(ax_h, esh, esp, res);
axis(ax_h, 'square')
ax_h.XLim = [min(esh(:)) max(esh(:))];
ax_h.YLim = [min(esp(:)) max(esp(:))];
set(ax_h, 'YDir', 'Normal')
cmap = bluered(65);
cmap(34:end, :) = [];
colormap(cmap)

plot(ax_h, esh, esp, 'color', 'k', 'linestyle', '-')
plot(ax_h, esh, 1.2*esh, 'color', 'k', 'linestyle', '--')



