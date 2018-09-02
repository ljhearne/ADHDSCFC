function [MAT,max_sz,max_szNull,nbsP] = NBS_StrucFuncBehav(behav,ADHDSC,CTRLSC,AllFC_AC,thresh,SCthresh,perms,Path,LOAD)
%explanation goes here

rng('shuffle')

Node = size(ADHDSC,1);
%% Analysis in real data

% perform SC-FC correlation for every node
[Fmatrix,~] = GLM_SCFCBEHAV(ADHDSC,CTRLSC,AllFC_AC,SCthresh,0,behav);
disp('Analysis - real data done');

%% Analysis - permutations
if LOAD == 0
    Fmatrix_null = zeros(Node,Node,perms);
    
    parfor p = 1:perms

        [FmatrixNull(:,:,p),~] = GLM_SCFCBEHAV(ADHDSC,CTRLSC,AllFC_AC,SCthresh,1,behav);
p
    end
    
    %save the data
    out = [Path,'Results/NBSSCFCBEHAV_SCthresh',num2str(SCthresh),'_Perms',num2str(perms),',.mat'];
    save(out,'FmatrixNull');
    
elseif LOAD == 1
    try
        load([Path,'Results/NBSSCFCBEHAV_SCthresh',num2str(SCthresh),'_Perms',num2str(perms),',.mat']);
    catch
        disp('NBS permutations file doesnt seem to exist');
    end
end

disp('Analysis - perms done');

%% Threshold and compare
% Code largely from the Network Based Statistic by Andrew Zalesky.
[MAT,max_sz] = NBS_findcomponents(Fmatrix,thresh);


for p = 1:perms
    [~,max_szNull(p)] = NBS_findcomponents(FmatrixNull(:,:,p),thresh);
end

%calculate significance (not sure if this is how NBS does it).
perc = prctile(max_szNull,1:0.001:100);
[~, p] = min(abs(perc'-max_sz));
nbsP = (length(perc)-p)/length(perc);
disp(['NBS p value = ',num2str(nbsP)]);


% % 
% % % plot component p -value
% % figure
% % hist(max_sz_null)
% % disp(['real largest component = ',num2str(max_sz)])
% % 
% % % plot original r-values
% % idx = triu(NBS_MAT,1);
% % 
% % figure
% % plot(r1(idx));hold on
% % plot(r2(idx),'k');hold on
% % 
% % % correlation plot (for paper)
% % [idx_row,idx_col] = find(sparse(idx));
% % 
% % for i = 1:length(idx_row)
% %     
% %     i_idx = idx_row(i);
% %     j_idx = idx_col(i);
% %     
% %     % ADHD
% %     x(:,i) = squeeze(AllFC_AC(i_idx,j_idx,1:N(2)));
% %     y(:,i) = squeeze(ADHDSC(i_idx,j_idx,:));
% %     
% %     %CTRL data
% %     x2(:,i) = squeeze(AllFC_AC(i_idx,j_idx,N(2)+1:end));
% %     y2(:,i) = squeeze(CTRLSC(i_idx,j_idx,:));
% %     
% % end
% % 
% % x_mean = mean(x,2);
% % y_mean = sum(y,2); % sum not mean
% % x2_mean = mean(x2,2);
% % y2_mean = sum(y2,2); % sum not mean
% % 
% % figure('Color',[1 1 1],'pos',[1000 600 400 225]);
% % figparam.lw = 1;
% % figparam.alpha = 0.5;
% % figparam.colctrl = [0.5 0.5 0.5];
% % figparam.coladhd = [1 0.5 0.1];
% % figparam.s = 30;
% % %-----------------------------------------------
% % subplot(1,2,1)
% % corr(x2_mean,y2_mean)
% % scatter(x2_mean,y2_mean,...
% %     'MarkerEdgeColor',figparam.colctrl,...
% %     'LineWidth',figparam.lw,...
% %     'MarkerFaceAlpha',figparam.alpha,...
% %     'MarkerEdgeAlpha',figparam.alpha,...
% %     'MarkerFaceColor', figparam.colctrl,...
% %     'SizeData',figparam.s);
% % hold on;
% % h = lsline;
% % set(h,'color',figparam.colctrl,'LineWidth',2)
% % 
% % xlabel('Mean FC');
% % ylabel('Mean SC');
% % set(gca,'FontName', 'Helvetica','FontSize', 10,'Box','off',...
% %     'TickDir','out','ygrid','off','XLim',[0.25 1],'YLim',[0 250]);
% % 
% % subplot(1,2,2)
% % corr(x_mean,y_mean)
% % scatter(x_mean,y_mean,...
% %     'MarkerEdgeColor',figparam.coladhd,...
% %     'LineWidth',figparam.lw,...
% %     'MarkerFaceAlpha',figparam.alpha,...
% %     'MarkerEdgeAlpha',figparam.alpha,...
% %     'MarkerFaceColor', figparam.coladhd,...
% %     'SizeData',figparam.s);
% % hold on;
% % h = lsline;
% % set(h,'color',figparam.coladhd,'LineWidth',2)
% % 
% % set(gca,'FontName', 'Helvetica','FontSize', 10,'Box','off',...
% %     'TickDir','out','ygrid','off','XLim',[0.25 1],'YLim',[0 250]);
% % 
% % % % anatomical vis
% % % figure
% % % load([basedir,'Luke/',ATLAS,'214Info/214parcellation_Yeo8Index.mat']);
% % % Comm = Yeo8Index; % network affiliations
% % %
% % % MNI = COG;
% % % mat2brainnet(NBS_MAT,MNI,Comm,ones(length(MNI),1));
% % % BrainNet_MapCfg([path2BN,'/Data/SurfTemplate/BrainMesh_Ch2.nv'],...
% % %     'NODE.node','EDGE.edge','BN_Hub_options.mat');
% % 
% % % save NBS results
% % outfile = ['NBSresult_perms',num2str(perms),...
% %     '_Zthresh',num2str(thresh),...
% %     '_SCthresh',num2str(SC_thresh),...
% %     '.mat'];
% % 
% % save(['Results/',outfile],'nbsP','z_matrix','z_matrix_thresh','Zthresh',...
% %     'perms','SC_thresh','NBS_MAT');
% % end