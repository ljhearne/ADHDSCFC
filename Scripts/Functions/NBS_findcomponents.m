function [MAT,max_sz] = NBS_findcomponents(Fmatrix,thresh)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% Code largely from the Network Based Statistic by Andrew Zalesky
% (NBSstats.m)
Node = size(Fmatrix,1);

%threshold the matrix
FmatrixThr = Fmatrix > thresh;

%get components within thresholded matrix
[comps,comp_sizes] = get_components(FmatrixThr); %BCT function

% size of components (greater than 1 & therefore useful)
ind_sz=find(comp_sizes>1);

sz_edges=zeros(1,length(ind_sz));
max_sz=0;
for i=1:length(ind_sz) %for each non trivial component
    
    n = find(ind_sz(i)==comps); % number of nodes implicated
    sz_edges(i) = sum(sum(FmatrixThr(n,n))); %count the number of edges (i.e., extent)
    
    if max_sz<sz_edges(i)
        max_sz = sz_edges(i); %max_sz equals new largest extent
        
        % save the edges to a Node x Node matrix for later
        idx = zeros(Node,Node);
        idx(n,n) = idx(n,n)+1;
        MAT = (idx+FmatrixThr)>1;
    end
end
end

