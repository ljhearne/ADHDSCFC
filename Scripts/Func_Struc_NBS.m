%explanation goes here

clearvars -except ADHDSC CTRLSC AllFC_AC
close all
rng('shuffle')
LASTN = maxNumCompThreads(40); % for server, apparently.

%--- inputs ------------%
LOAD = 1; % load the null perms from previous analysis
perms = 5000;
Zthresh = -1.96; %e.g., 1.96 = p = 0.05 ADHD > CTRL, 2.58
SC_thresh = 62; %how many SC data points are needed to consider a SC-FC correlation? 62 = 80% of smallest group

%--- paths -------------%
if ispc == 1
    basedir = 'L:/Lab_LucaC/4_ADHD_Luke_HY/';
else
    basedir = '/home/lucac/ldrive/Lab_LucaC/4_ADHD_Luke_HY/';
end

addpath([basedir,'Luke/x_scripts']); %path to custom scripts
addpath([basedir,'BCT/']); %path to BCT
addpath(genpath([basedir,'NBS1.2'])); %path to NBS
path2BN = [basedir,'BrainNet'];
addpath(genpath(path2BN)); % path to brainnet

ATLAS = 'Shaefer214/';
load([basedir,'Luke/',ATLAS,'214Info/Schaefer214_coordinates.mat']);

%--- load data ---------%
if exist('AllFC_AC','var')==1
    %do nothing
else
    load([basedir,'Luke/',ATLAS,'FC/AllFC_ADHD_CTRL.mat']);
    load([basedir,'Luke/',ATLAS,'SC/CTRLSC118.mat']);
    load([basedir,'Luke/',ATLAS,'SC/ADHDSC78.mat']);
end

SampSize(1) = size(CTRLSC,3); %sample size
SampSize(2) = size(ADHDSC,3);
Node = size(ADHDSC,1);
disp('data loaded');

%% Analysis - real data

% perform SC-FC correlation for every node
for i = 1:Node
    for j = 1:Node
        %ADHD data
        x = squeeze(AllFC_AC(i,j,1:SampSize(2)));
        y = squeeze(ADHDSC(i,j,:));
        
        %CTRL data
        x2 = squeeze(AllFC_AC(i,j,SampSize(2)+1:end));
        y2 = squeeze(CTRLSC(i,j,:));
        
        % exclude comparisons without enough data in the SC-matrices.
        if sum(y>0)<SC_thresh||sum(y2>0)<SC_thresh
            %assign 0
            r1(i,j) = 0;
            r2(i,j) = 0;
        else
            %remove any individuals with no SC data
            idx = y~=0;
            x = x(idx);
            y = y(idx);
            
            idx = y2~=0;
            x2 = x2(idx);
            y2 = y2(idx);
            
            %do correlation
            [r1(i,j),~] = corr(x,y);
            [r2(i,j),~] = corr(x2,y2);
            
            %test for difference - assign to matrix
            [z_matrix(i,j),~] = compare_coef(r1(i,j),r2(i,j),SampSize(2),SampSize(1));
        end
    end
end
disp('Analysis - real data done');

%% Analysis - permutations

if LOAD == 0
    z_matrix_null = zeros(Node,Node,perms);
    
    parfor p = 1:perms
        tic
        
        for i = 1:Node
            for j = 1:Node
                %ADHD data
                x = squeeze(AllFC_AC(i,j,1:SampSize(2)));
                y = squeeze(ADHDSC(i,j,:));
                
                %CTRL data
                x2 = squeeze(AllFC_AC(i,j,SampSize(2)+1:end));
                y2 = squeeze(CTRLSC(i,j,:));
                
                % permute data by shuffling group affiliation
                xp = [x;x2];
                yp = [y;y2];
                
                %shuffle index
                idx = randperm(length(xp));
                xp = xp(idx);
                yp = yp(idx);
                
                %redefine x's and y's
                x = xp(1:SampSize(2));
                y = yp(1:SampSize(2));
                x2 = xp(SampSize(2)+1:end);
                y2 = yp(SampSize(2)+1:end);
                
                % as before, exclude comparisons without enough data in the SC-matrices.
                if sum(y>0)<SC_thresh||sum(y2>0)<SC_thresh
                    %assign 0
                else
                    %remove any individuals with no SC data
                    idx = y~=0;
                    x = x(idx);
                    y = y(idx);
                    
                    idx = y2~=0;
                    x2 = x2(idx);
                    y2 = y2(idx);
                    
                    
                    %do correlation
                    [r1_null,~] = corr(x,y);
                    [r2_null,~] = corr(x2,y2);
                    % do test
                    [z_matrix_null(i,j,p),~] = compare_coef(r1_null,r2_null,SampSize(2),SampSize(1));
                end
            end
        end
        toc
    end
    
    outfile = ['z_matrix_null_isServer',num2str(isunix),'_perms',num2str(perms),...
        '_SCthresh',num2str(SC_thresh),'.mat'];
    save(outfile,'z_matrix_null');
    
elseif LOAD == 1
    outfile = ['z_matrix_null_isServer1','_perms',num2str(perms),...
        '_SCthresh',num2str(SC_thresh),'.mat'];
    load(outfile)
    
end

disp('Analysis - perms done');
%% Threshold and compare
% Code largely from NBS

%threshold the matrix
if Zthresh > 0
    z_matrix_thresh = z_matrix > Zthresh;
elseif Zthresh < 0
    z_matrix_thresh = z_matrix < Zthresh;
end

%get components from NBS
[comps,comp_sizes] = get_components(z_matrix_thresh);

ind_sz=find(comp_sizes>1); % size of components (greater than 1 & therefore useful)
sz_edges=zeros(1,length(ind_sz));

max_sz=0;
for i=1:length(ind_sz)
    n = find(ind_sz(i)==comps); % number of nodes implicated
    sz_edges(i) = sum(sum(z_matrix_thresh(n,n)))/2; %count the number of edges (i.e., extent)
    
    if max_sz<sz_edges(i)
        max_sz = sz_edges(i); %max_sz equals new largest extent
        
        % save the edges to a N x N matrix for later
        idx = zeros(Node,Node);
        idx(n,n) = idx(n,n)+1;
        NBS_MAT = (idx+z_matrix_thresh)>1;
    end
end

for p = 1:perms
    tmp = z_matrix_null(:,:,p);
    
    if Zthresh > 0
        tmp = tmp > Zthresh;
    elseif Zthresh < 0
        tmp = tmp < Zthresh;
    end
    
    [comps,comp_sizes] = get_components(tmp);
    XX(p) = max(comp_sizes);
    ind_sz=find(comp_sizes>1); % size of components (greater than 1 & therefore useful)
    sz_edges=zeros(1,length(ind_sz));
    
    max_sz_null(p)=0;
    for i=1:length(ind_sz)
        n = find(ind_sz(i)==comps); % number of nodes implicated
        sz_edges(i) = sum(sum(tmp(n,n)))/2; %count the number of edges (i.e., extent)
        
        if max_sz_null(p)<sz_edges(i)
            max_sz_null(p) = sz_edges(i); %max_sz equals new largest extent
        end
    end
    clear tmp
end

% calculate significance (do later, clearly significant)
perc = prctile(max_sz_null,1:0.001:100);
[~, p] = min(abs(perc'-max_sz));
nbsP = (length(perc)-p)/length(perc);
disp(['NBS p value = ',num2str(nbsP)]);
%% visualize component
close all

% plot component p -value
figure
hist(max_sz_null)
disp(['real largest component = ',num2str(max_sz)])

% plot original r-values
idx = triu(NBS_MAT,1);

figure
plot(r1(idx));hold on
plot(r2(idx),'k');hold on

% correlation plot (for paper)
[idx_row,idx_col] = find(sparse(idx));

for i = 1:length(idx_row)
    
    i_idx = idx_row(i);
    j_idx = idx_col(i);
    
    % ADHD
    x(:,i) = squeeze(AllFC_AC(i_idx,j_idx,1:SampSize(2)));
    y(:,i) = squeeze(ADHDSC(i_idx,j_idx,:));
    
    %CTRL data
    x2(:,i) = squeeze(AllFC_AC(i_idx,j_idx,SampSize(2)+1:end));
    y2(:,i) = squeeze(CTRLSC(i_idx,j_idx,:));
    
end

x_mean = mean(x,2);
y_mean = sum(y,2); % sum not mean
x2_mean = mean(x2,2);
y2_mean = sum(y2,2); % sum not mean

figure('Color',[1 1 1],'pos',[1000 600 400 225]);
figparam.lw = 1;
figparam.alpha = 0.5;
figparam.colctrl = [0.5 0.5 0.5];
figparam.coladhd = [1 0.5 0.1];
figparam.s = 30;
%-----------------------------------------------
subplot(1,2,1)
corr(x2_mean,y2_mean)
scatter(x2_mean,y2_mean,...
    'MarkerEdgeColor',figparam.colctrl,...
    'LineWidth',figparam.lw,...
    'MarkerFaceAlpha',figparam.alpha,...
    'MarkerEdgeAlpha',figparam.alpha,...
    'MarkerFaceColor', figparam.colctrl,...
    'SizeData',figparam.s);
hold on;
h = lsline;
set(h,'color',figparam.colctrl,'LineWidth',2)

xlabel('Mean FC');
ylabel('Mean SC');
set(gca,'FontName', 'Helvetica','FontSize', 10,'Box','off',...
    'TickDir','out','ygrid','off','XLim',[0.25 1],'YLim',[0 250]);

subplot(1,2,2)
corr(x_mean,y_mean)
scatter(x_mean,y_mean,...
    'MarkerEdgeColor',figparam.coladhd,...
    'LineWidth',figparam.lw,...
    'MarkerFaceAlpha',figparam.alpha,...
    'MarkerEdgeAlpha',figparam.alpha,...
    'MarkerFaceColor', figparam.coladhd,...
    'SizeData',figparam.s);
hold on;
h = lsline;
set(h,'color',figparam.coladhd,'LineWidth',2)

set(gca,'FontName', 'Helvetica','FontSize', 10,'Box','off',...
    'TickDir','out','ygrid','off','XLim',[0.25 1],'YLim',[0 250]);

% % anatomical vis
% figure
% load([basedir,'Luke/',ATLAS,'214Info/214parcellation_Yeo8Index.mat']);
% Comm = Yeo8Index; % network affiliations
% 
% MNI = COG;
% mat2brainnet(NBS_MAT,MNI,Comm,ones(length(MNI),1));
% BrainNet_MapCfg([path2BN,'/Data/SurfTemplate/BrainMesh_Ch2.nv'],...
%     'NODE.node','EDGE.edge','BN_Hub_options.mat');

% save NBS results
outfile = ['NBSresult_perms',num2str(perms),...
        '_Zthresh',num2str(Zthresh),...
        '_SCthresh',num2str(SC_thresh),...
        '.mat'];
    
save(['Results/',outfile],'nbsP','z_matrix','z_matrix_thresh','Zthresh',...
    'perms','SC_thresh','NBS_MAT');