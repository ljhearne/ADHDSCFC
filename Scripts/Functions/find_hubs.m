function [ConCount,ConStren,ConStrenNorm,MAT_out] = find_hubs(data,K)

%Code from Andrew Z. finds hub nodes and assigns labels to all connections
%in the matrix
% enter data for one participant (2d node by node matrix) and K-level for
% defining hub. Matrcies should be fully weighted and undirected (i.e. /2).

Nodes = size(data,1); %num of nodes in parc.
Nconn = (Nodes*(Nodes-1))/2;

d = sum(data); %weighted degree
%bet = betweenness_wei
ind = (d>=K); % binary list of hubs satisfying K
ind_ind = find(ind); % find hubs

%list hub-hubs
[u,v] = find(triu(data(ind,ind),1));
ind_rich = sub2ind([Nodes,Nodes],ind_ind(u),ind_ind(v));
ConCount(1) = length(ind_rich)/2; % number of connections

% list periphery (using the opposite index).
ind=~ind; ind_ind=find(ind);
[u,v]=find(triu(data(ind,ind),1));
ind_peri=sub2ind([Nodes,Nodes],ind_ind(u),ind_ind(v));
ConCount(3) = length(ind_peri)/2;

% calculate feeders (by subtraction)
ConCount(2) = Nconn - ConCount(3)- ConCount(2);

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

% Connectivity strength
ConStren(1) = sum(data(logical(MAT.hub)))/2; % total strength of connections
ConStren(2) = sum(data(logical(MAT.feed)))/2; % total strength of connections
ConStren(3) = sum(data(logical(MAT.peri)))/2; % total strength of connections

ConStrenNorm = ConStren ./ConCount;
end