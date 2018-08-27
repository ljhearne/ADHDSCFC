function draw_hubconnectome(MAT,COG,hubs,col,s,lw,profile)
% draws a basic top-down brain plot - MAT = hub connections in a 3d matrix,
% COG = MNI coords, hubs = index of hubs, col = colour,s=size, lw =
% linewidth,profile = 1'top' or 2'side'

%only take upper triangle
% MAT2 = zeros(size(MAT));
% MAT2 = triu(MAT,1);
% MAT2 = MAT2+MAT2';
% deg = sum(MAT2)/2;

% % draw top ~<thresh_lim> edges
% thresh = sort(abs(MAT2(:)),'descend');
% thresh = thresh(thresh_lim);
%
MAT2 = triu(MAT(:,:,1),1); % hubs only
c= 1;
for i = 1:length(MAT2)
    for j = 1:length(MAT2)
        
        if MAT2(i,j) ==1
            edge_x(c,:) = [COG(i,1),COG(j,1)];
            edge_y(c,:) = [COG(i,2),COG(j,2)];
            edge_z(c,:) = [COG(i,3),COG(j,3)];
            c=c+1;
        end
    end
end

if profile == 1
    
    % draw edges
    patchline(edge_x,edge_y,'edgecolor',col,'linewidth',2,'edgealpha',0.5); hold on
    clear edge_x edge_y
    
    % draw other nodes
    scatter(COG(:,1),COG(:,2),'SizeData',s,...
        'MarkerEdgeColor',col,...
        'MarkerFaceColor','w',...
        'MarkerEdgeAlpha', 0.5,...
        'MarkerFaceAlpha', 0.5,...
        'LineWidth',lw); hold on
    % draw hubs
    scatter(COG(hubs,1),COG(hubs,2),'SizeData',s*4,...
        'MarkerEdgeColor',col,...
        'MarkerFaceColor','w',...
        'MarkerEdgeAlpha', 1,...
        'MarkerFaceAlpha', 1,...
        'LineWidth',lw*3); hold on
    
elseif profile ==2
    
    patchline(edge_y,edge_z,'edgecolor',col,'linewidth',2,'edgealpha',0.5); hold on
    clear edge_x edge_y
    
    % draw other nodes
    scatter(COG(:,2),COG(:,3),'SizeData',s,...
        'MarkerEdgeColor',col,...
        'MarkerFaceColor','w',...
        'MarkerEdgeAlpha', 0.5,...
        'MarkerFaceAlpha', 0.5,...
        'LineWidth',lw); hold on
    % draw hubs
    scatter(COG(hubs,2),COG(hubs,3),'SizeData',s*4,...
        'MarkerEdgeColor',col,...
        'MarkerFaceColor','w',...
        'MarkerEdgeAlpha', 1,...
        'MarkerFaceAlpha', 1,...
        'LineWidth',lw*3); hold on
end
end

