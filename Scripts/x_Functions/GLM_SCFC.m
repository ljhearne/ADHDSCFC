function [stat_matrix,p_matrix] = GLM_SCFC(ADHDSC,CTRLSC,AllFC_AC,SCthresh,null)
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
    
    %ADHD data
    x = squeeze(AllFC_AC(i,j,1:N(2)));
    y = squeeze(ADHDSC(i,j,:));
    
    %CTRL data
    x2 = squeeze(AllFC_AC(i,j,N(2)+1:end));
    y2 = squeeze(CTRLSC(i,j,:));
    
    if null==1
        % this code shuffles the data for null permutations.
        
        % permute data by shuffling group affiliation
        xp = [x;x2];
        yp = [y;y2];
        
        %shuffle index
        idx = randperm(length(xp));
        xp = xp(idx);
        yp = yp(idx);
        
        %redefine x's and y's
        x = xp(1:N(2));
        y = yp(1:N(2));
        x2 = xp(N(2)+1:end);
        y2 = yp(N(2)+1:end);
        
    end
    
    % exclude comparisons without enough data in the SC-matrices.
    if sum(y>0)<SCthresh||sum(y2>0)<SCthresh
        
        %assign 0 and do nothing else
        stat_matrix(i,j) = 0;
        p_matrix(i,j) = 0;
    else
        
        %remove any individuals with no SC data
        idx = y~=0;
        x = x(idx);
        y = y(idx);
        
        idx = y2~=0;
        x2 = x2(idx);
        y2 = y2(idx);
        
        %test for difference GLM
        FC = [x;x2];
        SC = [y;y2];
        GROUP = [ones(length(x),1)*-1;ones(length(x2),1)];
        INT = GROUP .* SC;
        
        LM = fitlm([SC,GROUP,INT],FC,'linear');
        AN = anova(LM);
        stat_matrix(i,j) = AN{3,4};
        p_matrix(i,j) = AN{3,5};
    end
end
end

