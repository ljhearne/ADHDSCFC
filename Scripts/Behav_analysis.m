function behav = Behav_analysis(r,behav)
% HY's analysis (includes all subjects rather than just ADHD).

% All_Symp(:,1)~Inattention SNAP-IV (parent-rated);
%All_Symp(:,2)~Hyperactivity/Impulsivity SNAP-IV (parent-rated);
%All_Symp(:,3)~Inattention ASRS (self-rated);
%All_Symp(:,4)~Hyperactivity/Impulsivity ASRS (self-rated);

% All_Symp(1:78,:)~ADHD; All_Symp(79:196,:)~CTRL;
%behav = load([DataPath,'BrainBehavior.mat'],'All_Symp');
disp('---BEHAV statistics---');
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
%data = behav.All_Symp(idx,:);
data = behav.All_Symp_norm;
[behav.PCA.coeff,behav.PCA.score,behav.PCA.latent,behav.PCA.tsquared,behav.PCA.explained] = ...
    pca(data);
disp(['1st PCA factor accounts for ',num2str(behav.PCA.explained(1)),' % variance']);

% correlate with imaging (all SC-FC)
tmp = [r.all.ADHD,r.all.CTRL]';
[behav.PCACorr.r,behav.PCACorr.p] = corr(behav.PCA.score(:,1),tmp(idx), 'Type','Pearson');
disp(['Correlation with all SC-FC r = ',num2str(behav.PCACorr.r),', p = ',num2str(behav.PCACorr.p)])

[behav.r,behav.p]=corr(data,tmp(idx),'Type','Spearman'); %supplementary analysis (what if we don't do a PCA?)

% correlate with imaging (feeder only SC-FC)
tmp = [r.hub.ADHD(:,2);r.hub.CTRL(:,2)];
[behav.PCACorr.rhub,behav.PCACorr.phub] = corr(behav.PCA.score(:,1),tmp(idx), 'Type','Pearson');
disp(['Correlation with Feeder SC-FC r = ',num2str(behav.PCACorr.rhub),', p = ',num2str(behav.PCACorr.phub)])

[behav.rhub,behav.phub]=corr(data,tmp(idx),'Type','Spearman'); %supplementary analysis
disp(['Supplementary result: without PCA each individual behaviour is still',...
    ' negatively correlated with the SC-FC measure, r = ',...
    num2str(behav.rhub(1)),' ',num2str(behav.rhub(2)),' ',num2str(behav.rhub(3)),' ',num2str(behav.rhub(4)),...
    ' and p = ',...
    num2str(behav.phub(1)),' ',num2str(behav.phub(2)),' ',num2str(behav.phub(3)),' ',num2str(behav.phub(4)),...
    ]);
end

