function [sigma_noise, dist_var] = introduce_block_heterogeneity(node_attribute, mu_a, mu_b, std_a, std_b)
%%  Introduces variability in two types of nodes, the ones with the node attribute, 
%   and the ones without. 
%
% ARGUMENTS:
%    node_attribute    -- a vector of 0 and indicating which nodes have a certain (eg, hubs, nonhubs)  
%    mu_a              -- mean of the noise amplitude for nodes with node attribute
%    mu_b              -- mean of the noise amplitude for nodes witouth node attribute
%    var_a             -- std of the noise amplitude for nodes with node attribute
%    var_b             -- mean of the noise amplitude for nodes witouth node attribute
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
   
%}

% Get the diagonal of the sigma matrix
xx = zeros(length(node_attribute), 1);

% Node attribute is list of 0's and 1's
% The nodes that don't have that attribute, are set to the mean value mu_b
idx_ones  = find(node_attribute);
idx_zeros = find(~node_attribute);

% Make sure they are all positive
xx(idx_ones)  = sqrt((mu_a + std_a * randn(length(idx_ones), 1)).^2);
xx(idx_zeros) = sqrt((mu_b + std_b * randn(length(idx_zeros), 1)).^2);

% Recalculate the variance of all sigma_i
dist_var = var(xx);

% Return the modified sigma_noise matrix with sigma_i
sigma_noise = diag(xx);

end % function introduce_node_based_variability()