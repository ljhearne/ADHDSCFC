function [newdata] = normal_transform(data)

new_data = zeros(size(data));
for i = 1:size(data,2)
    
    rank = tiedrank(data(:,i));
    p = rank / ( length(rank) + 1 ); %# +1 to avoid Inf for the max point
    newdata(:,i) = norminv(p, 0, 1);
    
end

%We used a rank-based inverse Gaussian transformation29, to enforce
%Gaussianity for each of the SMs, producing S2. This transformation was
%used to avoid undue influence of potential outlier values, although we
%later confirmed that this normalisation made almost no difference to the
%final CCA results (see below).