function [stat_matrix,p_matrix] = GLM_SCFCBEHAV(ADHDSC,CTRLSC,AllFC_AC,SCthresh,null,behav)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N(1) = size(CTRLSC,3); %sample size
N(2) = size(ADHDSC,3);

Node = size(CTRLSC,1);
stat_matrix = zeros(Node,Node);
p_matrix = zeros(Node,Node);
idx = triu(ones(Node,Node),1);
[rowidx,colidx] = find(idx);

for edge = 1:length(rowidx)
    i = rowidx(edge);
    j = colidx(edge);
    
    %ADHD data only
    x = squeeze(AllFC_AC(i,j,1:N(2)));
    y = squeeze(ADHDSC(i,j,:));
    
    if null==1
        % this code shuffles the data for null permutations.
        
        % permute data by shuffling group affiliation
        xp = x;
        yp = y;
        
        %shuffle index
        idx = randperm(length(xp));
        xp = xp(idx);
        yp = yp(idx);
        
        %redefine x's and y's
        x = xp(1:N(2));
        y = yp(1:N(2));
        
    end
    
    % exclude comparisons without enough data in the SC-matrices.
    if sum(y>0)<SCthresh
        
        %assign 0 and do nothing else
        stat_matrix(i,j) = NaN;
        p_matrix(i,j) = NaN;
    else
        
        %remove any individuals with no SC data
        idx = y~=0;
        
        %test for difference GLM
        FC = x(idx);
        SC = y(idx);
        BEHAV1 = behav.Hypr(idx);
        BEHAV2 = behav.Inat(idx);
        INT = SC .* BEHAV1 .* BEHAV2;
        
        LM = fitlm([SC,INT],FC,'linear');
        AN = anova(LM);
        stat_matrix(i,j) = AN{2,4};
        p_matrix(i,j) = AN{2,5};
    end
end