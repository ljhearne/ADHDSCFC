function [deg,conCount,conStren,hubMat] = Struc_analysis(ADHDSC,CTRLSC,K)
% [deg,conCount,conStren,hubMat] = Struc_analysis(ADHDSC,CTRLSC,K)
% this function completes the main structural analysis for the ADHDSCFC
% paper. It calculates degree (weighted and binary for both groups) as well
% as connection class count (i.e., the sum of weights within each
% connection class) connection class strength (the count divided by the
% number of connections). It calculates hubs according to K (a percentage).
% The hub code is adapted from Andrew Zalesky.

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
disp(['Deg t-test, pval = ',num2str(P),' z = ',num2str(STATS.zval)]);
[P,~,STATS] = ranksum(deg.CTRLw,deg.ADHDw);
disp(['Weighted deg t-test, pval = ',num2str(P),' z = ',num2str(STATS.zval)]);


%% Hubs

%Klevel(1) = M+(1.5*S); %1.5 STD above mean


for i = 1:N(1)
    tmp = sort(sum(CTRLSC(:,:,i)),'descend');
    Klevel = tmp(round(length(tmp)*K));
    [conCount.CTRL(i,:),~,conStren.CTRL(i,:),hubMat.CTRL(:,:,:,i)] = find_hubs(CTRLSC(:,:,i),Klevel);
end

for i = 1:N(2)
    tmp = sort(sum(ADHDSC(:,:,i)),'descend');
    Klevel = tmp(round(length(tmp)*K));
    [conCount.ADHD(i,:),~,conStren.ADHD(i,:),hubMat.ADHD(:,:,:,i)] = find_hubs(ADHDSC(:,:,i),Klevel);
end

% Stats
disp('---Connection class statistics---');
for conn = 1:3
  %  [P,~,STATS] = ranksum(conCount.CTRL(:,conn),conCount.ADHD(:,conn));
  %  disp(['Hub COUNT    t-test- k = ',num2str(K),', conn = ',num2str(conn),' pval = ',num2str(P),' z = ',num2str(STATS.zval)]);
    
    [P,~,STATS] = ranksum(conStren.CTRL(:,conn),conStren.ADHD(:,conn));
    disp(['Hub STRENGTH t-test- k = ',num2str(K),', conn = ',num2str(conn),' pval = ',num2str(P),' z = ',num2str(STATS.zval)]);
end