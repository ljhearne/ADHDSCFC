clearvars
close all
clc

% This the main script for the ADHD structure-function paper. It is linked
% to the github repo: https://github.com/ljhearne/ADHDSCFC

% Please note:
% - This script assumes the data is organized in a very specific way which
% is outlined in the github repo
% - The basic NBS analyses are not incorporated in this script (e.g., healthy
% control versus ADHD structural matrices). I just used the GUI.
% - It is a little unclear what the "bandwidth" property does in the
% raincloud plots - may be worth checking.

%---------------------------------%
%%% Edit this %%%
%---------------------------------%

%Paths, functions and toolboxes
DataPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Data/';
DocsPath = '/Users/luke/Documents/Projects/ADHDStrucFunc/Docs/';

% Hub definition: top 15% of connections
K = 0.15;
addpath(genpath('Functions'));
addpath(genpath('Toolbox'));
Atlas = '214';

%---------------------------------%
%%% Stop editing %%%
%---------------------------------%

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
[~,r] = StrucFunc_analysis(ADHDSC,CTRLSC,AllFC_AC,hubMat);

%% Behaviour
% HY's analysis (includes all subjects rather than just ADHD).
% All_Symp(:,1)~Inattention SNAP-IV (parent-rated); 
%All_Symp(:,2)~Hyperactivity/Impulsivity SNAP-IV (parent-rated); 
%All_Symp(:,3)~Inattention ASRS (self-rated); 
%All_Symp(:,4)~Hyperactivity/Impulsivity ASRS (self-rated); 
  
% All_Symp(1:78,:)~ADHD; All_Symp(79:196,:)~CTRL;
behav = load([DataPath,'BrainBehavior.mat'],'All_Symp');

% '999' ~ no available data; 4th CTRL (All_Symp(82,:) ); 
idx = [1:81,83:196]; % exclude missing data

% are the distrubutions normal? (no)
for i = 1:4
    h=lillietest(behav.All_Symp(idx,i));
    if h == 1
        disp(['Distribution of values in behavioural variable ',num2str(i), ' is non-normal']);
    end
end

% transform the behavioural data for PCA
behav.All_Symp_log=log([behav.All_Symp(idx,1),behav.All_Symp(idx,2),...
    behav.All_Symp(idx,3),behav.All_Symp(idx,4)]);
behav.All_Symp_log(behav.All_Symp_log==-Inf)=0; % Inf values = 0

    %check that this transform worked (it doesn't);
    for i = 1:4
        h=lillietest(behav.All_Symp_log(:,i));
        if h == 1
            disp(['Log transformed values in var ',num2str(i), ' are still non-normal...']);
        end
    end

%try a different transform
for i = 1:4
    behav.All_Symp_norm(:,i) = normal_transform(behav.All_Symp(idx,i));
end

    %check that this transform worked (1/4)
    for i = 1:4
        h=lillietest(behav.All_Symp_norm(:,i));
        if h == 1
            disp(['Normal transformed values in var ',num2str(i), ' are still non-normal...']);
        end
    end

% Reduce dimesionality via PCA
data = behav.All_Symp(idx,:);
[behav.PCA.coeff,behav.PCA.score,behav.PCA.latent,behav.PCA.tsquared,behav.PCA.explained] = ...
    pca(data);
disp(['1st PCA factor accounts for ',num2str(behav.PCA.explained(1)),' % variance']);

% correlate with imaging (all SC-FC)
tmp = [r.all.ADHD,r.all.CTRL]';
[behav.PCACorr.r,behav.PCACorr.p] = corr(behav.PCA.score(:,1),tmp(idx), 'Type','Spearman');
disp(['Correlation with all SC-FC r = ',num2str(behav.PCACorr.r),', p = ',num2str(behav.PCACorr.p)])

[behav.r,behav.p]=corr(data,tmp(idx),'Type','Spearman'); %supplementary analysis (what if we don't do a PCA?)

% correlate with imaging (hubs only SC-FC)
tmp = [r.hub.ADHD(:,1);r.hub.CTRL(:,1)];
[behav.PCACorr.rhub,behav.PCACorr.phub] = corr(behav.PCA.score(:,1),tmp(idx), 'Type','Spearman');
disp(['Correlation with HUB SC-FC r = ',num2str(behav.PCACorr.rhub),', p = ',num2str(behav.PCACorr.phub)])

[behav.rhub,behav.phub]=corr(data,tmp(idx),'Type','Spearman'); %supplementary analysis 
disp(['Supplementary result: without PCA each individual behaviour is still',...
    ' negatively correlated with the SC-FC measure, r = ',...
    num2str(behav.rhub(1)),' ',num2str(behav.rhub(2)),' ',num2str(behav.rhub(3)),' ',num2str(behav.rhub(4)),...
    ' and p = ',...
    num2str(behav.phub(1)),' ',num2str(behav.phub(2)),' ',num2str(behav.phub(3)),' ',num2str(behav.phub(4)),...
    ]);

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

%draw raw structural matrices
close all
figure('Color','w','Position',[50 850 500 145]); hold on
setCmap(cols);

subplot(1,3,1)
data1 = mean(CTRLSC,3);
data1(logical(eye(size(data1)))) = 0;
i = prctile(data1(:),90);
imagesc(data1,[i*-1,i])
title('Control')
set(gca,'FontName', 'Helvetica','FontSize', 12);

subplot(1,3,2)
data = mean(ADHDSC,3);
data(logical(eye(size(data)))) = 0;
%data = data(idx,idx);
imagesc(data,[i*-1,i])
title('ADHD')
set(gca,'FontName', 'Helvetica','FontSize', 12);

subplot(1,3,3)
data = data1-data;
data(logical(eye(size(data)))) = 0;
%data = data(idx,idx);
imagesc(data,[i*-1,i])
title('Control - ADHD')
set(gca,'FontName', 'Helvetica','FontSize', 12);
saveas(gcf,[resultsdir,'Figure1a_matrices.svg']);


figure('Color','w','Position',[50 450 350 200]); hold on
% connectome density
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
figure('Color','w','Position',[450 450 600 200]); hold on
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
saveas(gcf,[resultsdir,'Figure2a_SChubclasses.svg']);
%% Figure 2b: Group hub topology
%figure('Color','w','Position',[850 450 450 450]); hold on

% have to calculate hubs at the group level
%CTRL
tmp = sort(sum(mean(CTRLSC,3)),'descend');
Klevel = tmp(round(length(tmp)*K));
tmp = sum(mean(CTRLSC,3));
ind = (tmp>=Klevel); % binary list of hubs satisfying K
ind = find(ind); % find hubs
[~,~,~,MAT] = find_hubs(mean(CTRLSC,3),Klevel);

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
tmp = sort(sum(mean(ADHDSC,3)),'descend');
Klevel = tmp(round(length(tmp)*K));
tmp = sum(mean(ADHDSC,3));
ind = (tmp>=Klevel); % binary list of hubs satisfying K
ind = find(ind); % find hubs
[~,~,~,MAT] = find_hubs(mean(ADHDSC,3),Klevel);

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
figure('Color','w','Position',[50 50 800 200]); hold on
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
saveas(gcf,[resultsdir,'Figure3_SCFChubclasses.svg']);
%% Figure 4: Association between SC-FC and behaviour (PCA)
figure('Color','w','Position',[900 25 300 250]); hold on

x = [r.hub.ADHD(:,1);r.hub.CTRL(:,1)];
x(82) = []; %remember! missing data
y = behav.PCA.score(:,1);

%ADHD only
scatter(x(1:78),y(1:78),...
    'MarkerEdgeColor',cl(2,:),...
    'MarkerFaceColor',cl(2,:),...
    'MarkerEdgeAlpha', 0.5,...
    'MarkerFaceAlpha', 0.5,...
    'LineWidth',1); hold on;
scatter(x(79:end),y(79:end),...
    'MarkerEdgeColor',cl(1,:),...
    'MarkerFaceColor',cl(1,:),...
    'MarkerEdgeAlpha', 0.5,...
    'MarkerFaceAlpha', 0.5,...
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