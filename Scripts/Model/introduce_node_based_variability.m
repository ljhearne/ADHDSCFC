function [sigma_noise, dist_std] = introduce_node_based_variability(sigma_noise, node_attribute, mu)
%%  Introduces variability in the nodes indicated in the index list
%
% ARGUMENTS:
%    sigma_noise       -- a diagonal matrix with values of sigma noise per node
%    node_attribute    -- a vector of 0 and indicating which nodes have a certain (eg, hubs, nonhubs)  
%    mu                -- mean of the noise amplitude 
% OUTPUT:
%    sigma_noise -- modified diagonal matrixdiagonal matrix 
%    dist_std    -- the new standard deviation of the sampled sigma_i distribution 
%
% AUTHOR:
% 
%     PSL 2019, QIMR Berghofer
%
% USAGE:
%{  
   
    S = get_sigma_noise(number_of_nodes, dist_std, seed, mu);
%}




if nargin < 3
    mu = 1;
end

% Get the diagonal of the sigma matrix
xx = diag(sigma_noise);
% Node attribute is list of 0's and 1's
% The nodes that don't have that attribute, are set to the mean value
xx(node_attribute == 0) = mu;

% Recalculate the standard deviation of the sigm
dist_std = std(xx);

% Return the modified sigma_noise matrix
sigma_noise = diag(xx);

end % function introduce_node_based_variability()