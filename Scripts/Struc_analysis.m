function [deg,conCount,conStren,hubMat,hublist] = Struc_analysis(ADHDSC,CTRLSC,K)
% this function completes the main structural analysis for the ADHDSCFC
% paper. It calculates degree (weighted and binary for both groups) as well
% as connection class 'count' (i.e., the sum of weights within each
% connection class) connection class strength (the count divided by the
% number of connections). It then performs a non-parametric t-test
% (ranksum.m) on the data.

% INPUTS:
% - ADHDSC: node x node x participant ADHD 3d structural matrices
% - CTRLSC: node x node x participant ADHD 3d structural matrices
% - K: The hub definition (e.g., 0.15)

N(1) = size(CTRLSC,3); %sample size
N(2) = size(ADHDSC,3);

%NConn = (size(ADHDSC,1)*(size(ADHDSC,1)-1));

%% Degree
for i = 1:N(1)
    t = CTRLSC(:,:,i);
    deg.CTRL(i) = sum(sum((t>0)))/2; %binary
    deg.CTRLw(i) = sum(sum(t))/2; %weighted
end

for i = 1:N(2)
    t = ADHDSC(:,:,i);
    deg.ADHD(i) = sum(sum((t>0)))/2;
    deg.ADHDw(i) = sum(sum(t))/2;
end

%stats
disp('---Degree statistics---');
[P,~,STATS] = ranksum(deg.CTRL,deg.ADHD);
STATS.r = STATS.zval / sqrt(sum(N));
disp(['Deg t-test, pval = ',num2str(P),' z = ',num2str(STATS.zval),' r = ',num2str(STATS.r)]);
[P,~,STATS] = ranksum(deg.CTRLw,deg.ADHDw);
STATS.r = STATS.zval / sqrt(sum(N));
disp(['Weighted deg t-test, pval = ',num2str(P),' z = ',num2str(STATS.zval),' r = ',num2str(STATS.r)]);


%% Hubs

%Klevel(1) = M+(1.5*S); %1.5 STD above mean
for i = 1:N(1)
    [conCount.CTRL(i,:),~,conStren.CTRL(i,:),hubMat.CTRL(:,:,:,i),hublist.CTRL(i,:)] = find_hubs(CTRLSC(:,:,i),K);
end

for i = 1:N(2)
    [conCount.ADHD(i,:),~,conStren.ADHD(i,:),hubMat.ADHD(:,:,:,i),hublist.ADHD(i,:)] = find_hubs(ADHDSC(:,:,i),K);
end

% Stats
disp('---Connection class statistics---');
for conn = 1:3
  %  [P,~,STATS] = ranksum(conCount.CTRL(:,conn),conCount.ADHD(:,conn));
  %  disp(['Hub COUNT    t-test- k = ',num2str(K),', conn = ',num2str(conn),' pval = ',num2str(P),' z = ',num2str(STATS.zval)]);
    [P,~,STATS] = ranksum(conStren.CTRL(:,conn),conStren.ADHD(:,conn));
    STATS.r = STATS.zval / sqrt(sum(N));
    disp(['Hub STRENGTH t-test- k = ',num2str(K),', conn = ',num2str(conn),' pval = ',num2str(P),' z = ',num2str(STATS.zval),' r = ',num2str(STATS.r)]);
end