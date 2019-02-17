function [node_var, node_mean] = guesstimate_ou_var_mean(data)
%% Returns the correlation coefficient between SC and FC
%
% ARGUMENTS:
%    data         -- an T x Nodes x subjects array with functional timeseries.
%                    By default the last time points is used to calculate
%                    the estimates because under the assumption that at
%                    T=end is the asymptotic limit of the process.
% OUTPUT:
%    ou_estimates -- node-wise estimates of the expected value and variance of an OU process. 
%                    The subjects are the realisations. 
%
% AUTHOR:
% 
%     PSL 2019, QIMR Berghofer
%
% USAGE:
%{
    
%}

% Asymptotic node-wise variance 
node_var  = var(squeeze(data(end, :, :)), [], 2);
% Asymptotic node-wise mean
node_mean = mean(squeeze(data(end, :, :)), 2); 
    
end % guesstimate_ou_var_mean()