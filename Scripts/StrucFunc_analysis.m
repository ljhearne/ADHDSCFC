function [r,rnorm] = StrucFunc_analysis(ADHDSC,CTRLSC,AllFC_AC,hubMat)
% [r,rlog] = StrucFunc_analysis(ADHDSC,CTRLSC,AllFC_AC,hubMat) this
% function calculate SC-FC correlations within individual for the whole
% connectome, then hubs/feeders/peripheries. Every non-0 SC value is added
% to a vector and correlated with a vector of FC within the same edges. So
% for each analysis a distribution of r-values is formed for each group which can be
% compared with a between-groups t-test.

% INPUTS:
% - ADHDSC: node x node x participant ADHD 3d structural matrices
% - CTRLSC: node x node x participant ADHD 3d structural matrices
% - AllFC_AC: All subjects FC matrices. This was entered as one matrix due
% to how the data was organised and no other reason
% - hubMat: a node x node x connection class x individual 4d binary matrix
% (see find_hubs.m)

N(1) = size(CTRLSC,3); %sample size
N(2) = size(ADHDSC,3);

% connectome-wide
for i = 1:N(1)
    SC = CTRLSC(:,:,i);
    FC = AllFC_AC(:,:,i+N(2));
    
    idx = tril(ones(size(SC)));
    SC(logical(idx)) = 0; %make lower tri = 0
    
    idx = SC(:) ~=0; %index not zero values
    SCn = SC(idx);
    FCn = FC(idx);
    
    [r.all.CTRL(i)] = corr(SCn,FCn);
    %[rlog.all.CTRL(i)] = corr(log(SCn),FCn);
    [rnorm.all.CTRL(i)] = corr(normal_transform(SCn),FCn);
end

for i = 1:N(2)
    SC = ADHDSC(:,:,i);
    FC = AllFC_AC(:,:,i);
    
    idx = tril(ones(size(SC)));
    SC(logical(idx)) = 0; %make lower tri = 0
    
    idx = SC(:) ~=0; %index not zero values
    SCn = SC(idx);
    FCn = FC(idx);
    
    [r.all.ADHD(i)] = corr(SCn,FCn);
    %[rlog.all.ADHD(i)] = corr(log(SCn),FCn);
    [rnorm.all.ADHD(i)] = corr(normal_transform(SCn),FCn);
end

disp('---STRUC-FUNC statistics---');
%[P,~,STATS] = ranksum(r.all.CTRL,r.all.ADHD);
%disp(['Connectome-wide t-test, pval = ',num2str(P),' z = ',num2str(STATS.zval)]);
%[P,~,STATS] = ranksum(rlog.all.CTRL,rlog.all.ADHD);
%disp(['Connectome-wide LOG t-test, pval = ',num2str(P),' z = ',num2str(STATS.zval)]);
[P,~,STATS] = ranksum(rnorm.all.CTRL,rnorm.all.ADHD);
STATS.r = STATS.zval / sqrt(sum(N));
disp(['Connectome-wide t-test (NORMAL), pval = ',num2str(P),' z = ',num2str(STATS.zval),' r = ',num2str(STATS.r)]);
normality_tests(rnorm.all.CTRL,rnorm.all.ADHD);

% same analysis but divided by connection class.
for connType = 1:3
    for i = 1:N(1)
        SC = CTRLSC(:,:,i);
        FC = AllFC_AC(:,:,i+N(2));
        
        idx = tril(ones(size(SC)));
        SC(logical(idx)) = 0; %make lower tri = 0
        
        % index only the appropriate connection type
        idx = ~hubMat.CTRL(:,:,connType,i);
        SC(logical(idx)) = 0; %make connections not of interest = 0
        
        idx = SC(:) ~=0; %index not zero values
        SCn = SC(idx);
        FCn = FC(idx);
        
        [r.hub.CTRL(i,connType)] = corr(SCn,FCn);
        %[rlog.hub.CTRL(i,connType)] = corr(log(SCn),FCn);
        [rnorm.hub.CTRL(i,connType)] = corr(normal_transform(SCn),FCn);
    end
    
    for i = 1:N(2)
        SC = ADHDSC(:,:,i);
        FC = AllFC_AC(:,:,i);
        
        idx = tril(ones(size(SC)));
        SC(logical(idx)) = 0; %make lower tri = 0
        
        % index only the appropriate connection type
        idx = ~hubMat.ADHD(:,:,connType,i);
        SC(logical(idx)) = 0; %make connections not of interest = 0
        
        idx = SC(:) ~=0; %index not zero values
        SCn = SC(idx);
        FCn = FC(idx);
        
        [r.hub.ADHD(i,connType)] = corr(SCn,FCn);
        %[rlog.hub.ADHD(i,connType)] = corr(log(SCn),FCn);
        [rnorm.hub.ADHD(i,connType)] = corr(normal_transform(SCn),FCn);
    end
end

classlabel = {'Hub','Feeder','Periphery'};
for connType = 1:3
    %[P,~,STATS] = ranksum(r.hub.CTRL(:,connType),r.hub.ADHD(:,connType));
    %disp([classlabel{connType},' t-test, pval = ',num2str(P),' z = ',num2str(STATS.zval)]);
    %[P,~,STATS] = ranksum(rlog.hub.CTRL(:,connType),rlog.hub.ADHD(:,connType));
    %disp([classlabel{connType},' LOG t-test, pval = ',num2str(P),' z = ',num2str(STATS.zval)]);
    [P,~,STATS] = ranksum(rnorm.hub.CTRL(:,connType),rnorm.hub.ADHD(:,connType));
    STATS.r = STATS.zval / sqrt(sum(N));
    disp([classlabel{connType},' NORMAL t-test, pval = ',num2str(P),', pbonf(3) = ',num2str(P*3),', z = ',num2str(STATS.zval),' r = ',num2str(STATS.r)]);
    normality_tests(rnorm.hub.CTRL(:,connType),rnorm.hub.ADHD(:,connType));
end
end