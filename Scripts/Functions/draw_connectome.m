function draw_connectome(MAT,COG,s,thresh_lim)
%draws a very basic brain plot of nodes/edges:draw_connectome(MAT,COG,s,thresh_lim)

%only take upper triangle
MAT2 = zeros(size(MAT));
MAT2 = triu(MAT,1);
MAT2 = MAT2+MAT2';
deg = sum(MAT2)/2;

% draw top ~<thresh_lim> edges
thresh = sort(abs(MAT2(:)),'descend');
thresh = thresh(thresh_lim);

c= 1;
for i = 1:length(MAT2)
    for j = 1:length(MAT2)
        if abs(MAT2(i,j)) >= thresh
            if MAT2(i,j) > 0
                edge_x_pos(c,:) = [COG(i,1),COG(j,1)];
                edge_y_pos(c,:) = [COG(i,2),COG(j,2)];
                edge_z_pos(c,:) = [COG(i,3),COG(j,3)];
                c=c+1;
            elseif MAT2(i,j) < 0
                edge_x_neg(c,:) = [COG(i,1),COG(j,1)];
                edge_y_neg(c,:) = [COG(i,2),COG(j,2)];
                edge_z_neg(c,:) = [COG(i,3),COG(j,3)];
                c=c+1;
            end
        end
    end
end

%draw edges
try
patchline(edge_x_pos,edge_y_pos,'edgecolor','r','linewidth',1,'edgealpha',0.2); hold on
end
try
patchline(edge_x_neg,edge_y_neg,'edgecolor','b','linewidth',1,'edgealpha',0.2); hold on
end
% draw nodes
scatter(COG(:,1),COG(:,2),[],deg,'SizeData',s); hold on

end

