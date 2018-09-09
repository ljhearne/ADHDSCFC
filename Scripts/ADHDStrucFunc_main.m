clearvars
close all

% This the main script for the ADHD structure-function paper.


% TO DO:
% - Incorporate NBS into the script if possible (allows for better
%replicability, we already know the results).
% - Understand bandwidth property (affects visualization, not the actual
%results)

%Paths, functions and toolboxes
DataPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Data/';
DocsPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/';
addpath(genpath('x_Functions'));
addpath(genpath('Toolbox'));
%path2BN = [pwd,'/x_Toolbox/BrainNet'];

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

% Colours for plots
[cb] = cbrewer('qual','Set3',12,'pchip');
cl(1,:) = [0.5 0.5 0.5];
cl(2,:) = cb(4,:);

%% Structural analysis
K = 0.15; % Hub definition: top 15% of connections
[deg,conCount,conStren,hubMat] = Struc_analysis(ADHDSC,CTRLSC,K);

%% Structure - function analysis
% by connection class - hub, feeder & periphery.
[r,logr] = StrucFunc_analysis(ADHDSC,CTRLSC,AllFC_AC,hubMat);

% Thought to test by network but this might introduce methodological
% concerns regarding how confident we can be in correlations derived from
% smaller numbers.

%% Behaviour
behav.Inat = behav.raw(1:N(2),15);
behav.Hypr = behav.raw(1:N(2),16);
[behav.r,behav.p] = corr([behav.Inat,behav.Hypr,r.all.ADHD',r.hub.ADHD]);

disp('---BEHAV statistics---');
disp('no correlation between behaviour and struc-func');

%% Exploratory analysis SC-FC NBS
% this analysis fundamentally differs from previous as we are now doing the
% correlation/regression ACROSS subjects, rather than within subjects.
%SCthresh = round(N(2)*.80);
% SCthresh = 50;
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
Stab.perms = 1000;
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
figure('Color','w','Position',[50 450 350 350]); hold on

% add imagesc plots ordered by FC network.
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

%% Figure 2: Structural connection classes
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

%% Figure 3: Group hub topology
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

%% Figure 4: SC-FC correlations
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

line([.3 .3],[.25*-15 .75*-15],'Color','k');
text(.31,.5*-15,'*','HorizontalAlignment','center','FontSize',16)

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
    
    if i ==1
        line([.5 .5],[.25*-7 .75*-7],'Color','k');
        text(.51,.5*-7,'*','HorizontalAlignment','center','FontSize',16)
    end
end

%% Figure 5: no association between SC-FC and behaviour
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

%% Figure 6: Stability test
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