function sigma_noise = get_sigma_noise(number_of_nodes, distribution_std, seed, mu)
%% Returns a diagonal matrix with the values of sigma (standard deviation)
%  of the noise. The  function also checks if the values of sigma are negative, 
%  which they shouldn't be as they represent standard deviation.   
%
% ARGUMENTS:
%    number_of_nodes       -- scalar with the number of nodes in the SC or FC matrix.
%    distribution_variance -- scalar with the standard deviation of the distribution (Normal)  
%                             used to generate the nonuniform sigma values
% OUTPUT:
%    sigma_noise -- diagonal matrix 
%
% AUTHOR:
% 
%     PSL 2019, QIMR Berghofer
%
% USAGE:
%{  
    seed     = 33;
    dist_std = 2^-1;
    number_of_nodes = 42;
    mu = 1;
    S = get_sigma_noise(number_of_nodes, dist_std, seed, mu);
%}

% NOTE: This function does not check if the values of sigma are negative, 
%        which they shouldn't be as they represent standard deviation.       


if nargin < 3
    seed = 42;
    mu = 1;
end

if nargin < 4
    mu = 1;
end

rng(seed)

% Strategy for setting nonuniform sigma based on distribution
sigma_noise = diag(mu + distribution_std * randn(number_of_nodes, 1));

if any(sigma_noise < 0)
    disp(['There are negative values in the  matrix. Update mu or distribution_dist.' ...
         'Converting to a folded normal distribution'])
    sigma_noise = abs(sigma_noise);

end
end % function get_sigma_vec()