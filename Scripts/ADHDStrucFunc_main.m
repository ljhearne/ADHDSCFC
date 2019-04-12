clearvars
close all
clc

% This the main script for the ADHD structure-function paper. It is linked
% to the github repo: https://github.com/ljhearne/ADHDSCFC. For information
% about specific analyses see the functions.

% Please note:
% - This script assumes the data is organized in a very specific way which
% is outlined in the github repo
% - The basic NBS analyses are not incorporated in this script (e.g., healthy
% control versus ADHD structural matrices). I just used the GUI.

%---------------------------------%
%---------------------------------%

%Paths, functions and toolboxes
DataPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Data/';
DocsPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/';

% Hub definition: top 15% of connections
K = 0.15;
addpath(genpath('Functions'));
addpath(genpath('Toolbox'));

% Parcellation
Atlas = 'Schaefer214';
%Atlas = 'Shen268';
%Atlas = 'Brainnetome_246';
FIGS = 0; %draw figures?
SUPP = 0; %do supplemental analyses?

%---------------------------------%
%---------------------------------%
% Load data
load([DataPath,Atlas,'/',Atlas,'_Info/',Atlas,'_coordinates.mat']);
load([DataPath,Atlas,'/',Atlas,'_Info/',Atlas,'parcellation_Yeo8Index.mat']);
load([DataPath,Atlas,'/','SC/CTRLSC118.mat']);
load([DataPath,Atlas,'/','SC/ADHDSC78.mat']);
load([DataPath,Atlas,'/','FC/AllFC_ADHD_CTRL.mat']);

N(1) = size(CTRLSC,3); %sample size
N(2) = size(ADHDSC,3);

% Brainetome ROIs 117,118 should be deleted due to poor fMRI data in those
% regions
if strcmp(Atlas,'Brainnetome_246')
    tmp = [117,118];
    COG(tmp,:) =      [];
    Yeo8Index(tmp) =  [];
    
    ADHDSC(tmp,:,:) = [];
    ADHDSC(:,tmp,:) = [];
    
    CTRLSC(tmp,:,:) = [];
    CTRLSC(:,tmp,:) = [];
    
    AllFC_AC(tmp,:,:) = [];
    AllFC_AC(:,tmp,:) = [];
elseif strcmp(Atlas,'Shen268')
    COG(223:end,:)       =  [];
    Yeo8Index(223:end,:) =  [];
end

disp('Data is loaded');

%results directory
resultsdir = [DocsPath,'Results/K',num2str(K*100),'/',Atlas,'/'];
mkdir(resultsdir);
delete([resultsdir,'stats.txt']);
diary([resultsdir,'stats.txt']);

% Colours for plots
[cb] = cbrewer('qual','Set3',12,'pchip');
cl(1,:) = [0.5 0.5 0.5];
cl(2,:) = cb(1,:);

cols = [0, 0, 144
    255, 255, 255
    144, 0, 0]./255;

cols = [cl(2,:)
    1,1,1
    cl(1,:)];

%% Structural analysis
[deg,conCount,conStren,hubMat,hublist] = Struc_analysis(ADHDSC,CTRLSC,K);

%% Structure - function analysis
% by connection class - hub, feeder & local.
[~,r] = StrucFunc_analysis(ADHDSC,CTRLSC,AllFC_AC,hubMat);

%% Behaviour
% HY's analysis (includes all subjects)

%for reference
%All_Symp(:,1)~Inattention SNAP-IV (parent-rated);
%All_Symp(:,2)~Hyperactivity/Impulsivity SNAP-IV (parent-rated);
%All_Symp(:,3)~Inattention ASRS (self-rated);
%All_Symp(:,4)~Hyperactivity/Impulsivity ASRS (self-rated);

behav = load([DataPath,'BrainBehavior.mat'],'All_Symp');
behav = Behav_analysis(r,behav);

%% Functional follow up
% Feeder's seem to be the most interesting connection class. Here are are
% doing a follow up analysis where we split each feeder SC-FC vector into
% functional networks and redo the analysis.

feed = Func_followup(CTRLSC,ADHDSC,AllFC_AC,Yeo8Index,hubMat);

diary off

%% Supplementary Tables
if SUPP ==1
    
    % Supplementary Table 1
    Stable = zeros(length(COG),9);
    Stable(:,1) = 1:length(COG);
    Stable(:,2:4) = round(COG,1);
    Stable(:,5) = Yeo8Index;
    
    %hubs - group
    [~,~,~,~,ind] = find_hubs(mean(CTRLSC,3),K);
    Stable(:,6) = ind;
    
    [~,~,~,~,ind] = find_hubs(mean(ADHDSC,3),K);
    Stable(:,8) = ind;
    
    %hubs - individual
    ind = (sum(hublist.CTRL,1)/N(1))*100;
    Stable(:,7) = round(ind,1);
    
    ind = (sum(hublist.ADHD,1)/N(2))*100;
    Stable(:,9) = round(ind,1);
    
    l = {'Node','X' 'Y' 'Z' 'Network' 'GroupC' 'IndividualC' 'GroupA' 'IndividualA'};
    Stable = array2table(Stable,'VariableNames',l);
    writetable(Stable,[resultsdir,'SuppTable1.csv']);
    
    % Supplementary Table 2
    Stable = zeros(10,3);
    paths = {[DocsPath,'Results/K',num2str(K*100),'/Schaefer214/stats.txt'],...
        [DocsPath,'Results/K',num2str(K*100),'/Shen268/stats.txt'],...
        [DocsPath,'Results/K',num2str(K*100),'/Brainnetome_246/stats.txt']};
    for i = 1:length(paths)
        
        Stable(:,i) = Supplement_stats(paths{i},K);
    end
    
    l = {'Schaefer214','Shen268','Brainnetome246'};
    r = {'Degree (p)','Weighted degree (p)','SC Hubs (p)','SC Feeder (p)','SC Periphery (p)',...
        'SC-FC (p)','SC-FC Hubs (p)','SC-FC Feeders (p)','SC-FC Periphery (p)',...
        'Feeder-behaviour correlation (r)'};
    Stable = array2table(Stable,'VariableNames',l,'RowNames',r);
    writetable(Stable,[resultsdir,'SuppTable2.xlsx'],'WriteRowNames',1);
    
    % Supplementary Table 3
    Stable = zeros(10,3);
    paths = {[DocsPath,'Results/K12.5/',Atlas,'/stats.txt'],...
        [DocsPath,'Results/K15/',Atlas,'/stats.txt'],...
        [DocsPath,'Results/K17.5/',Atlas,'/stats.txt']};
    nK = [.125, .15,.175];
    for i = 1:length(paths)
        Stable(:,i) = Supplement_stats(paths{i},nK(i));
    end
    
    l = {'hub125','hub15','hub175'};
    Stable = array2table(Stable,'VariableNames',l,'RowNames',r);
    writetable(Stable,[resultsdir,'SuppTable3.xlsx'],'WriteRowNames',1);
    
    %Supplementary Table 4
    Stable = zeros(5,4);
    Stable(1:4,:) = round(behav.PCA.coeff,2);
    Stable(5,:) = round(behav.PCA.explained,2);
    l = {'Component1','Component2','Component3','Component4',};
    r = {'Inattention SNAP-IV (parent-rated)',...
        'Hyperactivity/Impulsivity SNAP-IV (parent-rated)',...
        'Inattention ASRS (self-rated)',...
        'Hyperactivity/Impulsivity ASRS (self-rated)',...
        'Variance explained'};
    
    Stable = array2table(Stable,'VariableNames',l,'RowNames',r);
    writetable(Stable,[resultsdir,'SuppTable4.xlsx'],'WriteRowNames',1);
    
    %Supplementary Table 5
    Stable = [];
    Stable = round(behav.rhub,2);
    Stable(:,2) = round(behav.phub,4);
    l = {'Correlation','p'};
    r(5) = [];
    
    Stable = array2table(Stable,'VariableNames',l,'RowNames',r);
    writetable(Stable,[resultsdir,'SuppTable5.xlsx'],'WriteRowNames',1);
    
    % Supplementary stability test
    Stab.perms = 100; % I would reccomend 1000
    Stab.subs = 50;
    
    for perms = 1:Stab.perms
        for i = 1:Stab.subs
            
            tmp1=r.hub.CTRL(:,2);
            tmp2=r.hub.ADHD(:,2);
            
            % delete 'i' subjects from each group
            idx = randperm(length(tmp1));
            tmp1(idx(1:i))=[];
            
            idx = randperm(length(tmp2));
            tmp2(idx(1:i))=[];
            
            [~,~,STATS] = ranksum(tmp1,tmp2);
            Stab.zval(i,perms) = STATS.zval;
        end
    end
end


%% Figure 1: Structural degree and weighted degree

if FIGS ==1
    
    %draw raw structural matrices
    close all
    figure('Color','w','Position',[50 850 400 240]); hold on
    setCmap(cols);
    
    subplot(2,3,1)
    data1 = mean(AllFC_AC(:,:,N(2)+1:end),3);
    data1(logical(eye(size(data1)))) = 0;
    i = 1.2;
    imagesc(data1,[i*-1,i])
    title('Control')
    set(gca,'FontName', 'Helvetica','FontSize', 12);
    
    subplot(2,3,2)
    data = mean(AllFC_AC(:,:,1:N(2)),3);
    data(logical(eye(size(data)))) = 0;
    %data = data(idx,idx);
    imagesc(data,[i*-1,i])
    title('ADHD')
    set(gca,'FontName', 'Helvetica','FontSize', 12);
    
    subplot(2,3,3)
    data = data1-data;
    data(logical(eye(size(data)))) = 0;
    %data = data(idx,idx);
    imagesc(data,[i*-1,i])
    title('Control - ADHD')
    set(gca,'FontName', 'Helvetica','FontSize', 12);
    %colorbar
    
    subplot(2,3,4)
    data1 = mean(CTRLSC,3);
    data1(logical(eye(size(data1)))) = 0;
    i = 10;
    imagesc(data1,[i*-1,i])
    %title('Control')
    set(gca,'FontName', 'Helvetica','FontSize', 12);
    
    subplot(2,3,5)
    data = mean(ADHDSC,3);
    data(logical(eye(size(data)))) = 0;
    %data = data(idx,idx);
    imagesc(data,[i*-1,i])
    %title('ADHD')
    set(gca,'FontName', 'Helvetica','FontSize', 12);
    
    subplot(2,3,6)
    data = data1-data;
    data(logical(eye(size(data)))) = 0;
    %data = data(idx,idx);
    imagesc(data,[i*-1,i])
    %title('Control - ADHD')
    set(gca,'FontName', 'Helvetica','FontSize', 12);
    %colorbar
    saveas(gcf,[resultsdir,'Figure1a_matrices.svg']);
    
    % connectome density
    figure('Color','w','Position',[50 450 350 200]); hold on
    
    subplot(1,2,1)
    title('')
    [~, ~, u] = ksdensity(deg.CTRL);
    h1 = raincloud_plot('X',deg.CTRL,'color', cl(1,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(1,:),...
        'box_dodge',1,'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1,...
        'bandwidth',u);
    h2 = raincloud_plot('X',deg.ADHD,'color', cl(2,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(2,:),...
        'box_dodge',1,'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1,...
        'bandwidth',u);
    legend([h1{1} h2{1}], {'Control', 'ADHD'},'Location','best')
    xlabel('Degree');
    set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off','view',[90 -90],'Ytick',[]);
    title('Density')
    %set(gca,'Xtick',0:3000:12000);
    hAxes = gca;
    hAxes.XAxis.Exponent = 2;
    subplot(1,2,2)
    title('')
    [~, ~, u] = ksdensity(deg.CTRLw);
    h1 = raincloud_plot('X',deg.CTRLw,'color', cl(1,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(1,:),...
        'box_dodge',1,'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1,...
        'bandwidth',u);
    h2 = raincloud_plot('X',deg.ADHDw,'color', cl(2,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(2,:),...
        'box_dodge',1,'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1,...
        'bandwidth',u);
    
    xlabel('Weighted degree');
    set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off','view',[90 -90],'Ytick',[]);
    %set(gca,'Xtick',0:5000:100000);
    hAxes = gca;
    hAxes.XAxis.Exponent = 3;
    saveas(gcf,[resultsdir,'Figure1b_degree.svg']);
    
    %% Figure 2a: Structural connection classes
    figure('Color','w','Position',[450 450 500 150]); hold on
    data = conStren; % choose whether to vis connectivity weighted/nonweighted
    
    classlabel = {'Hub','Feeder','local'};
    for i = 1:3
        subplot(1,3,i)
        title('')
        
        h1 = raincloud_plot('X',data.CTRL(:,i),'color', cl(1,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(1,:),...
            'box_dodge',1,'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1);
        h2 = raincloud_plot('X',data.ADHD(:,i),'color', cl(2,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(2,:),...
            'box_dodge',1,'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1);
        
        xlabel('Mean strength');
        ylabel(classlabel{i})
        yl = get(gca, 'YLim');
        set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off','view',[90 -90],'Ytick',[]);
        
        if i == 1
            set(gca,'Xtick',0:2:10);
            set(gca,'Xlim',[2,10]);
        elseif i ==2
            set(gca,'Xlim',[1,4]);
        elseif i ==3
            set(gca,'Xlim',[1,5]);
            set(gca,'YLim',[yl(1),yl(2)+.1]);
        end
    end
    saveas(gcf,[resultsdir,'Figure2a_SChubclasses.svg']);
    
    %% Hub topology - individual level
    
    figure('Color','w','Position',[50 850 300 50]); hold on
    
    [nets,idx] = sort(Yeo8Index,'ascend');
    a = sum(hublist.CTRL,1)/N(1); a = a(idx);
    b = sum(hublist.ADHD,1)/N(2); b = b(idx);
    data = [a',b'];
    colorbar
    i = 1;
    
    figure('Color','w','Position',[50 850 300 150]); hold on
    % actual data
    setCmap(cols);
    imagesc(data',[i*-1,i]);
    netIDS = find(diff(nets));
    netIDS(end+1) = length(nets);
    for i = 1:length(netIDS)
        h = line([netIDS(i) netIDS(i) ],[.25 3.5]);
        set(h,'Color','k');
    end
    set(gca,'FontName', 'Helvetica','FontSize', 12,'ylim',[.25 2.5],'xlim',[0 length(nets)]);
    set(gca,'YTick',1:2,'YTickLabel', {'Control' 'ADHD'});
    set(gca,'XTick',netIDS,'XTickLabel', {'Vis' 'SomMat' 'DorstAttn' 'SalVentAttn' 'Limbic' 'Control' 'Default' 'Misc'});
    xtickangle(45);
    saveas(gcf,[resultsdir,'Figure3c_indHubs.svg']);
    
    %% Hub topology - group level
    % have to calculate hubs at the group level
    %CTRL
    [~,~,~,MAT,ind] = find_hubs(mean(CTRLSC,3),K);
    
    MATn = zeros(size(MAT(:,:,1)));
    for i = 1:3
        idx = logical(MAT(:,:,i));
        MATn(idx) = 4-i;
    end
    
    net = ones(length(MAT),1);
    Nsize = net;
    Nsize(ind) = 5; %hubs
    out = [resultsdir,'Hub_CTRL'];
    
    surficeEdgeVis(MATn,COG,net,Nsize,[1,3],out)
    
    % % ADHD
    [~,~,~,MAT,ind] = find_hubs(mean(ADHDSC,3),K);
    
    MATn = zeros(size(MAT(:,:,1)));
    for i = 1:3
        idx = logical(MAT(:,:,i));
        MATn(idx) = 4-i;
    end
    
    net = ones(length(MAT),1);
    Nsize = net;
    Nsize(ind) = 5; %hubs
    out = [resultsdir,'Hub_ADHD'];
    
    surficeEdgeVis(MATn,COG,net,Nsize,[1,3],out)
    %% Figure 3: SC-FC correlations
    figure('Color','w','Position',[50 50 667 150]); hold on
    data = r; % choose 'logr' or 'r'
    xlims = [0 .55]; %same across plots for comparison
    
    subplot(1,4,1)
    title('')
    [~, ~, u] = ksdensity(data.all.CTRL);
    h1 = raincloud_plot('X',data.all.CTRL,'color', cl(1,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(1,:),...
        'box_dodge',1,'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1,...
        'bandwidth',u);
    h2 = raincloud_plot('X',data.all.ADHD,'color', cl(2,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(2,:),...
        'box_dodge',1,'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1,...
        'bandwidth',u);
    
    xlabel('r value');
    ylabel('Connectome');
    set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off','view',[90 -90],'Ytick',[]);
    
    set(gca,'Xtick',0:0.1:1);
    %set(gca,'YLim',[-15 15]); % Need to verify this code.
    set(gca,'Xlim',[0,0.5]);
    
    classlabel = {'Hub','Feeder','local'};
    for i = 1:3
        subplot(1,4,1+i)
        title('')
        [~, ~, u] = ksdensity(data.hub.CTRL(:,i));
        h1 = raincloud_plot('X',data.hub.CTRL(:,i),'color', cl(1,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(1,:),...
            'box_dodge',1,'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1);
        h2 = raincloud_plot('X',data.hub.ADHD(:,i),'color', cl(2,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(2,:),...
            'box_dodge',1,'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1);
        
        xlabel('r value');
        ylabel(classlabel{i});
        set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off','view',[90 -90],'Ytick',[]);
        set(gca,'Xlim',[0,0.5]);
        yl = get(gca, 'YLim');
        
        if i == 3
            set(gca,'YLim',[yl(1),yl(2)+1.7]);
        end
    end
    saveas(gcf,[resultsdir,'Figure3_SCFChubclasses.svg']);
    saveas(gcf,[resultsdir,'Figure3_SCFChubclasses.jpeg']);
    %% Figure 4: Association between SC-FC and behaviour (PCA)
    figure('Color','w','Position',[900 25 200 150]); hold on
    
    x = [r.hub.ADHD(:,2);r.hub.CTRL(:,2)];
    x(82) = []; %remember! missing data
    y = behav.PCA.score(:,1);
    size_data = 20;
    %ADHD only
    scatter(x(1:78),y(1:78),...
        'MarkerEdgeColor',cl(2,:),...
        'MarkerFaceColor',cl(2,:),...
        'MarkerEdgeAlpha', 0.7,...
        'MarkerFaceAlpha', 0.7,...
        'SizeData',size_data,...
        'LineWidth',1); hold on;
    scatter(x(79:end),y(79:end),...
        'MarkerEdgeColor',cl(1,:),...
        'MarkerFaceColor',cl(1,:),...
        'MarkerEdgeAlpha', 0.5,...
        'MarkerFaceAlpha', 0.5,...
        'SizeData',size_data,...
        'LineWidth',1); hold on;
    
    ax = scatter(x,y,...
        'MarkerEdgeColor','k',...
        'MarkerEdgeAlpha', 0,...
        'LineWidth',1); hold on;
    
    h = lsline;
    
    delete(h(2)); delete(h(3));
    set(h(1),'LineWidth',.1, 'Color','k');
    
    set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
    ylabel('PCA component loading');
    xlabel('SC-FC correlation');
    
    saveas(gcf,[resultsdir,'Figure4_behav.svg']);
    %     %% Figure 5: Stability test
    %     figure('Color','w','Position',[900 25 400 200]); hold on
    %
    %     for i = 1:Stab.subs
    %         d = Stab.zval(i,:);
    %         CI = ConfInt(d);
    %         %CI = [min(d),max(d)];
    %
    %         if min(CI) > 1.96
    %             line([i,i],[CI(1),CI(2)],'Color','k');
    %             scatter(i,mean(d),'k','filled');
    %         else
    %             line([i,i],[CI(1),CI(2)],'Color',[0.5 0.5 0.5]);
    %             scatter(i,mean(d),'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
    %         end
    %     end
    %
    %     % real data point
    %     [~,~,STATS] = ranksum(r.hub.CTRL(:,1),r.hub.ADHD(:,1));
    %     scatter(0,STATS.zval,'MarkerFaceColor',cl(2,:),'MarkerEdgeColor',cl(2,:));
    %
    %     set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
    %     xlim([0 50]);
    %     ylabel('Z-score');
    %     xlabel('Subjects deleted');
    %     saveas(gcf,[resultsdir,'Figure5_stabilityTest.jpeg']);
    %
    %% Feeder follow up
    figure('Color','w','Position',[50 850 225 160]); hold on
    
   % ncols = [1, 1, 1
   %     .85, .93, .98
   %     cl(2,1), cl(2,2), cl(2,3)];
    
    %setCmap(ncols);
    
    tmp = mean(feed.CTRL,3);
    tmp2 = mean(feed.ADHD,3);
    data = tmp - tmp2;
    imagesc(flipud(data),[0,0.026]);
    colorbar
    
    %title('Control > ADHD')
    %set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
    set(gca,'XTick',1:3,'XTickLabel', {'Control' 'Default-mode' 'Sensory'});
    xtickangle(45);
    set(gca,'YTick',1:3,'YTickLabel', {'Sensory' 'Default-mode' 'Control'});
    saveas(gcf,[resultsdir,'Feederfollowup.svg']);
    
    %% Methods figure
    figure('Color','w','Position',[50 850 120 270]); hold on
    setCmap(cols);
    
    [~,idx] = sort(Yeo8Index,'ascend'); %sort by FC networks
    sub = 1;
    subplot(2,1,1)
    data = CTRLSC(:,:,sub);
    data(logical(eye(size(data)))) = 0;
    i = prctile(data(:),90);
    imagesc(data,[i*-1,i])
    title('Sub01 Structure')
    set(gca,'FontName', 'Helvetica','FontSize', 12);
    set(gca,'Xtick',[]);
    set(gca,'Ytick',[]);
    
    subplot(2,1,2)
    data = AllFC_AC(:,:,sub);
    data(logical(eye(size(data)))) = 0;
    i = max(max(max(data)));
    imagesc(data,[i*-1,i])
    title('Sub01 Function')
    set(gca,'FontName', 'Helvetica','FontSize', 12);
    set(gca,'Xtick',[]);
    set(gca,'Ytick',[]);
    saveas(gcf,[resultsdir,'Methods_matrices.svg']);
    
    figure('Color','w','Position',[50 850 150 120]); hold on
    data = CTRLSC(:,:,sub);
    idx = tril(ones(size(data)),sub);
    data(logical(idx)) = 0;
    idx = data>0;
    y = normal_transform(data(idx));
    x = AllFC_AC(:,:,1);
    x = x(idx);
    
    scatter(x,y,'Marker','.','SizeData',15,...
        'MarkerEdgeColor',cl(1,:),'MarkerEdgeAlpha',0.5);
    ylabel('Structure');
    xlabel('Function');
    h = lsline;
    set(h,'LineWidth',2,'Color','k');
    saveas(gcf,[resultsdir,'Methods_correlation.svg']);
    
    %% Supplementary 1: results of the log transform.
    figure('Color','w','Position',[50 850 800 400]); hold on
    
    data = CTRLSC(:,:,1);
    idx = tril(ones(size(data)),sub);
    
    n = 5;
    for i = 1:n
        data = ADHDSC(:,:,i);
        data(logical(idx)) = 0;
        idx = data>0;
        y = data(idx);
        
        subplot(3,n,i)
        histogram(y)
        
        subplot(3,n,n+i)
        histogram(log(y))
        
        subplot(3,n,n*2+i)
        histogram(normal_transform(y))
    end
    saveas(gcf,[resultsdir,'Example_transform.jpeg']);
    
    %% Supplementary 2: results of transforming behaviour
    % note that this makes no difference to the results (improves relationship
    % between PCA and hubs very slightly)
    figure('Color','w','Position',[50 850 800 400]); hold on
    for i = 1:4
        subplot(3,4,i)
        hist(behav.All_Symp([1:81,83:end],i)); %missing! data
        title('non transformed');
        
        subplot(3,4,i+4)
        hist(behav.All_Symp_log(:,i));
        title('log transformed');
        
        subplot(3,4,i+8)
        hist(behav.All_Symp_norm(:,i));
        title('normal transformed');
    end
    saveas(gcf,[resultsdir,'Example_behaviour_transform.jpeg']);
end