function [ConCount,ConStren,ConStrenNorm,MAT_out,hublist] = find_hubs(data,K)

%This code idenitifies hubs and then the different connection classes, you
%input a node x node matrix and the K (number of nodes) level. It uses 4
%different graph metrics to define hubs (see below).

Nodes = size(data,1); %num of nodes in parc.
Nconn = (Nodes*(Nodes-1))/2;

%% Adaptation of Hsiang-Yuan's code to detect hubs across graph metrics

%degree
degree = sum(data>0);            %do graph metric
degree_rank = tiedrank(degree);  %rank graph metric
ranks(:,1) = degree_rank;        %concatenate to matrix

%weighted degree
weighted_degree(:,1) = sum(data);
weighted_degree_rank = tiedrank(weighted_degree);
ranks(:,2) = weighted_degree_rank; 

%betweeness
betweeness(:,1) = betweenness_wei(data);
betweeness_rank = tiedrank(betweeness);
ranks(:,3) = betweeness_rank; 

%nodal communicability
nodal_communicability(:,1) = subgraph_centrality(data);
nodal_communicability = transpose(nodal_communicability);
nodal_communicability_rank = tiedrank(nodal_communicability);
ranks(:,4) = nodal_communicability_rank;

%combine the metrics
AVGranks = mean(ranks,2); %average across the ranks
sortAVGranks = sortrows(AVGranks,'descend');

%% Adaptation of Andrew's code to organise hubs, feeders and periphery connections

%find and index the hubs (top K % of measure defined above
Klevel = sortAVGranks(floor(K*Nodes));
ind = (AVGranks>=Klevel);
ind_ind = find(ind);
hublist=ind; %save for later

%list hub-hubs
[u,v] = find(triu(data(ind,ind),1));
ind_rich = sub2ind([Nodes,Nodes],ind_ind(u),ind_ind(v));
%ConCount(1) = length(ind_rich)/2; % number of connections

% list periphery (using the opposite index).
ind=~ind; ind_ind=find(ind);
[u,v]=find(triu(data(ind,ind),1));
ind_peri=sub2ind([Nodes,Nodes],ind_ind(u),ind_ind(v));
%ConCount(3) = length(ind_peri)/2;

% calculate feeders (by subtraction)
%ConCount(2) = Nconn - ConCount(3)- ConCount(2);

% rich club calculation
MAT.hub = zeros(Nodes,Nodes);
MAT.peri = zeros(Nodes,Nodes);
MAT.feed = zeros(Nodes,Nodes);

MAT.hub(ind_rich)=1;
MAT.peri(ind_peri)=1;
MAT.feed = triu(ones(Nodes,Nodes),1)-MAT.hub-MAT.peri;

MAT_out = zeros(Nodes,Nodes,3);
MAT_out(:,:,1) = MAT.hub;
MAT_out(:,:,2) = MAT.feed;
MAT_out(:,:,3) = MAT.peri;

% Count the number of edges in each class accounting for sparsity of the
% matrix.
idx = data>0;
for i = 1:3
	tmp = idx(logical(MAT_out(:,:,i)));  
    ConCount(i) = sum(tmp);
end

% Connectivity strength (summed)
ConStren(1) = sum(data(logical(MAT.hub))); % total strength of connections
ConStren(2) = sum(data(logical(MAT.feed))); % total strength of connections
ConStren(3) = sum(data(logical(MAT.peri))); % total strength of connections

ConStrenNorm = ConStren ./ConCount;
end