function feed = Func_followup(CTRLSC,ADHDSC,AllFC_AC,net,hubMat)
% function to complete the follow-up feeder analysis. 
% - takes edges from specific networks (defined below) and correlates
% between structure and function data (like the rest of the paper). Only
% does this correlation if there are at least 50 data points, if not,
% returns a NAN value.

% To ensure we have enough edges in our networks we add together all the
% control networks.
new_net = zeros(size(net));
idx = net==3; new_net(idx) = 1; % Task-positive
idx = net==4; new_net(idx) = 1;
idx = net==6; new_net(idx) = 1;
idx = net==7; new_net(idx) = 2; % Default-mode
idx = net==1; new_net(idx) = 3; % Sensory
idx = net==2; new_net(idx) = 3;

N(1) = size(CTRLSC,3); %sample size
N(2) = size(ADHDSC,3);

% Controls
for i = 1:N(1)
    
    feedidx = hubMat.CTRL(:,:,2,i); %index of feeders for this p
    
    for j = 1:max(new_net)
        for k = 1:max(new_net)
            
            netidx = zeros(size(feedidx));
            
            net1 = new_net==j; %network index j
            net2 = new_net==k; %network index k
            
            netidx(net1,net2) = 1;
            netidx(net2,net1) = 1;
            
            % find the common between feeders and network idx
            new_idx = [];
            new_idx = netidx + feedidx == 2;
            
            SC = CTRLSC(:,:,i);
            SC = SC(new_idx);
            FC = AllFC_AC(:,:,i+N(2));
            FC = FC(new_idx);
            
            idx = SC(:) ~=0; %index non zero values
            SCn = SC(idx);
            FCn = FC(idx);
            
            if length(SCn) > 50
                [feed.CTRL(j,k,i)] = corr(normal_transform(SCn),FCn);
            else
                [feed.CTRL(j,k,i)] = NaN;
            end
        end
    end
end

% ADHD
for i = 1:N(2)
    
    feedidx = hubMat.ADHD(:,:,2,i); %index of feeders for this p
    
    for j = 1:max(new_net)
        for k = 1:max(new_net)
            netidx = zeros(size(feedidx));
            
            net1 = new_net==j;
            net2 = new_net==k;
            
            netidx(net1,net2) = 1;
            netidx(net2,net1) = 1;
            
            % find the common between feeders and network idx
            new_idx = [];
            new_idx = netidx + feedidx == 2;
            
            SC = ADHDSC(:,:,i);
            SC = SC(new_idx);
            FC = AllFC_AC(:,:,i);
            FC = FC(new_idx);
            
            idx = SC(:) ~=0; %index not zero values
            SCn = SC(idx);
            FCn = FC(idx);
            
            if length(SCn) > 50
                [feed.ADHD(j,k,i)] = corr(normal_transform(SCn),FCn);
            else
                [feed.ADHD(j,k,i)] = NaN;
            end
        end
    end
end

% test the NANS before doing this analysis
disp('Number of NaNs in feeder matrix across all participants:');
disp(sum(isnan(cat(3,feed.CTRL,feed.ADHD)),3));

% Do the statistical tests.
for j = 1:max(new_net)
    for k = 1:max(new_net)
        x = squeeze(feed.CTRL(j,k,:));
        y = squeeze(feed.ADHD(j,k,:));
        [feed.pmat(j,k),~,STATS] = ranksum(x,y);
        feed.zmat(j,k) = STATS.zval;
    end
end