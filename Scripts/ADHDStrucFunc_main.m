clearvars
close all
clc
% This the main script for the ADHD structure-function paper.


% Please note:
% The basic NBS analyses are not incorporated in this script (e.g., healthy
% control versus ADHD structural matrices).
% It is a little unclear what the "bandwidth" property does in the
% raincloud plots - may be worth checking.


%%% You need to edit this %%%
%Paths, functions and toolboxes
DataPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Data/';
DocsPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/';
% Hub definition: top 15% of connections
K = 0.15;
%%%%%%%% Stop editing %%%%%%%

addpath(genpath('x_Functions'));
addpath(genpath('Toolbox'));

% Atlas
Atlas = '214';

% Load data
load([DataPath,'Schaefer',Atlas,'/',Atlas,'Info/Schaefer',Atlas,'_coordinates.mat']);
load([DataPath,'Schaefer',Atlas,'/',Atlas,'Info/',Atlas,'parcellation_Yeo8Index.mat']);
load([DataPath,'Schaefer',Atlas,'/','SC/CTRLSC118.mat']);
load([DataPath,'/Schaefer',Atlas,'/','SC/ADHDSC78.mat']);
load([DataPath,'/Schaefer',Atlas,'/','FC/AllFC_ADHD_CTRL.mat']);
[behav.raw] = xlsread([DataPath,'/AdultADHD_FS_Cov_20180510.xlsx']);

N(1) = size(CTRLSC,3); %sample size
N(2) = size(ADHDSC,3);

%results directory
resultsdir = [DocsPath,'Results/K',num2str(K*100),'/',Atlas,'/'];
mkdir(resultsdir);
diary([resultsdir,'stats.txt']);

% Colours for plots
[cb] = cbrewer('qual','Set3',12,'pchip');
cl(1,:) = [0.5 0.5 0.5];
cl(2,:) = cb(4,:);

cols = [0, 0, 144
     255, 255, 255
     144, 0, 0]./255;

%% Structural analysis
[deg,conCount,conStren,hubMat] = Struc_analysis(ADHDSC,CTRLSC,K);

%% Structure - function analysis
% by connection class - hub, feeder & periphery.
r = StrucFunc_analysis(ADHDSC,CTRLSC,AllFC_AC,hubMat);

%% Behaviour
behav.Inat = behav.raw(1:N(2),15);
behav.Hypr = behav.raw(1:N(2),16);
[behav.r,behav.p] = corr([behav.Inat,behav.Hypr,r.all.ADHD',r.hub.ADHD]);

disp('---BEHAV statistics---');
disp(['Inat rval = ',num2str(behav.r(1,3:end)),' (all, hub, feed, peri)']);
%disp(['Inat rval = ',num2str(behav.p(1,3:end)),' (all, hub, feed, peri)']);

disp(['Hyper rval = ',num2str(behav.r(2,3:end)),' (all, hub, feed, peri)']);
%disp(['Hyper rval = ',num2str(behav.p(2,3:end)),' (all, hub, feed, peri)']);
disp('no consistent correlations');
%% Exploratory analysis SC-FC NBS (commented out)
% this analysis fundamentally differs from previous as we are now doing the
% correlation/regression ACROSS subjects, rather than within subjects.

%SCthresh = round(N(2)*.80);

% perms = 100;
% F = 7.5;
% %[MAT,max_sz,max_szNull,nbsP] = NBS_StrucFunc(ADHDSC,CTRLSC,AllFC_AC,F,...
% %    SCthresh,perms,DocsPath,1);
%
% SCthresh = 78;
% perms = 100;
% F = 5;
%
% [NBSMAT,max_sz,max_szNull,nbsP] = NBS_StrucFuncBehav(behav,ADHDSC,CTRLSC,AllFC_AC,F,...
%     SCthresh,perms,DocsPath,1);

%% Stability test
Stab.perms = 100; % I would reccomend 1000
Stab.subs = 50;
for perms = 1:Stab.perms
    for i = 1:Stab.subs
        
        tmp1=r.hub.CTRL(:,1);
        tmp2=r.hub.ADHD(:,1);
        
        % delete 'i' subjects from each group
        idx = randperm(length(tmp1));
        tmp1(idx(1:i))=[];
        
        idx = randperm(length(tmp2));
        tmp2(idx(1:i))=[];
        
        [~,~,STATS] = ranksum(tmp1,tmp2);
        Stab.zval(i,perms) = STATS.zval;
    end
end

%% Figure 1: Structural degree and weighted degree
diary off

close all
figure('Color','w','Position',[50 850 350 150]); hold on
setCmap(cols);

[~,idx] = sort(Yeo8Index,'ascend'); %sort by FC networks

subplot(1,2,1)
data = mean(AllFC_AC(:,:,1:N(2)),3);
data(logical(eye(size(data)))) = 0;
data = data(idx,idx);
imagesc(data,[max(max(max(data)))*-1,max(max(max(data)))])

subplot(1,2,2)
data = mean(AllFC_AC(:,:,N(2)+1:end),3);
data(logical(eye(size(data)))) = 0;
data = data(idx,idx);
imagesc(data,[max(max(max(data)))*-1,max(max(max(data)))])
saveas(gcf,[resultsdir,'Figure1a_matrices.jpeg']);

figure('Color','w','Position',[50 450 350 350]); hold on

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
xlabel('Summed degree');
set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off','view',[90 -90],'Ytick',[]);
set(gca,'Xtick',0:3000:12000);

subplot(1,2,2)
title('')
[~, ~, u] = ksdensity(deg.CTRLw);
h1 = raincloud_plot('X',deg.CTRLw,'color', cl(1,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(1,:),...
    'box_dodge',1,'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1,...
    'bandwidth',u);
h2 = raincloud_plot('X',deg.ADHDw,'color', cl(2,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(2,:),...
    'box_dodge',1,'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1,...
    'bandwidth',u);

xlabel('Summed weighted degree');
set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off','view',[90 -90],'Ytick',[]);
set(gca,'Xtick',0:5000:100000);
saveas(gcf,[resultsdir,'Figure1b_degree.jpeg']);
%% Figure 2a: Structural connection classes
figure('Color','w','Position',[450 450 525 350]); hold on
data = conStren; % choose whether to vis connectivity weighted/nonweighted

classlabel = {'hub','feeder','periphery'};
for i = 1:3
    subplot(1,3,i)
    title('')
    [~, ~, u] = ksdensity(data.CTRL(:,i));
    h1 = raincloud_plot('X',data.CTRL(:,i),'color', cl(1,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(1,:),...
        'box_dodge',1,'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1,...
        'bandwidth',u);
    h2 = raincloud_plot('X',data.ADHD(:,i),'color', cl(2,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(2,:),...
        'box_dodge',1,'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1,...
        'bandwidth',u);
    
    xlabel(['Mean ',classlabel{i},' strength']);
    set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off','view',[90 -90],'Ytick',[]);
    %set(gca,'Xtick',0:3000:12000);
end
saveas(gcf,[resultsdir,'Figure2a_SChubclasses.jpeg']);
%% Figure 2b: Group hub topology
figure('Color','w','Position',[850 450 525 525]); hold on

% have to calculate hubs at the group level
%CTRL
tmp = sort(sum(mean(CTRLSC,3)),'descend');
Klevel = tmp(round(length(tmp)*K));
tmp = sum(mean(CTRLSC,3));
ind = (tmp>=Klevel); % binary list of hubs satisfying K
ind = find(ind); % find hubs
[~,~,~,MAT] = find_hubs(mean(CTRLSC,3),Klevel);

subplot(3,2,[1 3])
draw_hubconnectome(MAT,COG,ind,cl(1,:),20,1,1);
axis off
title('CTRL Hub topology');

subplot(3,2,5)
draw_hubconnectome(MAT,COG,ind,cl(1,:),20,1,2);
axis off

% ADHD
tmp = sort(sum(mean(ADHDSC,3)),'descend');
Klevel = tmp(round(length(tmp)*K));
tmp = sum(mean(ADHDSC,3));
ind = (tmp>=Klevel); % binary list of hubs satisfying K
ind = find(ind); % find hubs
[~,~,~,MAT] = find_hubs(mean(ADHDSC,3),Klevel);

subplot(3,2,[2 4])
draw_hubconnectome(MAT,COG,ind,cl(2,:),20,1,1);
axis off
title('ADHD Hub topology');

subplot(3,2,6)
draw_hubconnectome(MAT,COG,ind,cl(2,:),20,1,2);
axis off
saveas(gcf,[resultsdir,'Figure2b_Hubs.jpeg']);
%% Figure 3: SC-FC correlations
figure('Color','w','Position',[50 50 800 350]); hold on
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

legend([h1{1} h2{1}], {'Control', 'ADHD'},'Location','best')
xlabel('r value');
ylabel('Connectome');
set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off','view',[90 -90],'Ytick',[]);

set(gca,'Xtick',0:0.1:1);
set(gca,'YLim',[-15 15]); % Need to verify this code.

classlabel = {'Hub','Feeder','Periphery'};
for i = 1:3
    subplot(1,4,1+i)
    title('')
    [~, ~, u] = ksdensity(data.hub.CTRL(:,i));
    h1 = raincloud_plot('X',data.hub.CTRL(:,i),'color', cl(1,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(1,:),...
        'box_dodge',1,'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'line_width',1,...
        'bandwidth',u);
    h2 = raincloud_plot('X',data.hub.ADHD(:,i),'color', cl(2,:),'box_on',1,'alpha',0.5,'cloud_edge_col', cl(2,:),...
        'box_dodge',1,'box_dodge_amount', .75, 'dot_dodge_amount', .75, 'box_col_match', 0,'line_width',1,...
        'bandwidth',u);
    
    xlabel('r value');
    ylabel(classlabel{i});
    set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off','view',[90 -90],'Ytick',[]);
    set(gca,'Xtick',0:0.1:1);
end
saveas(gcf,[resultsdir,'Figure3_SCFChubclasses.jpeg']);
%% Figure 4: no association between SC-FC and behaviour
figure('Color','w','Position',[900 25 400 200]); hold on

subplot(1,2,1)
scatter(r.all.ADHD',behav.Inat,...
    'MarkerEdgeColor','k',...
    'MarkerEdgeAlpha', 0.5,...
    'LineWidth',1); hold on;
h = lsline;
set(h,'LineWidth',1, 'Color',cl(2,:));
set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
ylabel('Inattention score');
xlabel('SC-FC correlation');

subplot(1,2,2)
scatter(r.all.ADHD',behav.Hypr,...
    'MarkerEdgeColor','k',...
    'MarkerEdgeAlpha', 0.5,...
    'LineWidth',1); hold on;
h = lsline;
set(h,'LineWidth',1, 'Color',cl(2,:));
set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
ylabel('Hyperactivity score');
xlabel('SC-FC correlation');
saveas(gcf,[resultsdir,'Figure4_behav.jpeg']);
%% Figure 5: Stability test
figure('Color','w','Position',[900 25 400 200]); hold on

for i = 1:Stab.subs
    d = Stab.zval(i,:);
    CI = ConfInt(d);
    %CI = [min(d),max(d)];
    
    if min(CI) > 1.96
        line([i,i],[CI(1),CI(2)],'Color','k');
        scatter(i,mean(d),'k','filled');
    else
        line([i,i],[CI(1),CI(2)],'Color',[0.5 0.5 0.5]);
        scatter(i,mean(d),'MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
    end
end

% real data point
[~,~,STATS] = ranksum(r.hub.CTRL(:,1),r.hub.ADHD(:,1));
scatter(0,STATS.zval,'MarkerFaceColor',cl(2,:),'MarkerEdgeColor',cl(2,:));

set(gca,'FontName', 'Helvetica','FontSize', 12,'box','off');
xlim([0 50]);
ylabel('Z-score');
xlabel('Subjects deleted');
saveas(gcf,[resultsdir,'Figure5_stabilityTest.jpeg']);