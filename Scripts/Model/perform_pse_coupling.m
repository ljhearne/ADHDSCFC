%   This function performs a parameter space exploration overvalues of coupling strength. 
%   for each value it calculates the analytic FC (PCC) from a structural connectivity matrix
%   
%   This function computes the analytical functional connectivity, measured
%   as Pearson Correlation Coefficient, starting from structural data.
%   The analytical operation is obtained under the assumption of an
%   underlying stochastic linear model. The parameter of the model, the global
%   coupling c, is set in order to optimize the fit between the analytical
%   and the empirical FC.
%
%   Inputs:
%   empSC -->   is the n x n structural connectivity matrix, where n is the
%               number of regions in the parcellation
%   empPCC -->  is the n x n functional connectivity matrix, measured as PCC
%   prec -->    the algorithm works until the critical point c_{critic}. The
%               parameter space is explored from c=0 to c=c_{critic}-prec.
%   visualize_figure --> type '1' if you want to see a figure of aSC
%
%   Output:
%   aPCC --> is a matrix n x n whose entries represent an estimate of the
%           functional connectivity for the value of the parameter which
%           maximizes the correlation with empirical data
%
% Original from Saggio et al 2016. 
% Modified by PSL 2019 QIMR to output more data and handle nonunfirom sigma

function [aPCC, c_best, c_critic, corr_max, c_range, corr_vec] = perform_pse_coupling(empSC, empPCC, pre_c, sigma_noise)

if nargin < 4
    sigma_noise= 1;
end

number_of_nodes = size(empSC, 1);  % number of regions in the parcellation

% define the parameter range to be explored
c_critic = 1/max(eig(empSC));  % find critical value for the parameter after
                               % which the system destabilizes (at the extact
                               % critical value the inverse matrix cannot be
                               % computed)
                               
cmax = c_critic - c_critic*pre_c;   % PSL: explore a bit below the crritical value 
cmin = 2^-9;                        % PSL: set a nonzero minimum

% discretization of parameter range
step_number = 1024;                 %PSL: NOTE: harcoded number
step_size   = (cmax-cmin)/step_number;

% explore parameter range
corr_vec = zeros(1,step_number);  % will save the correlation between empPCC and aPCC
                                  % for each value of the parameter global_coupling considered
c_range = zeros(step_number+1, 1);

for ii=1:(step_number+1)
    c_range(ii) = cmin+(ii-1)*step_size;
    aPCC_c      = calculate_analyticalFCPCC(empSC, c_range(ii), sigma_noise);  
    % Compute the correlation removing the diagonal elements. For some
    % reason the original used to include the diagonal.
    corr_vec(ii)  = corr(empPCC(tril(ones(number_of_nodes), -1) > 0), aPCC_c(tril(ones(number_of_nodes), -1) > 0));
end

% find best value for global_coupling
[~, index] = max(corr_vec);              % find highest correlation
c_best     = cmin+(index-1)*step_size;   % find best value for coupling in our range

% compute best aPCC, that is, the one that has the highest corrrelation with
% the empirical one
aPCC     = calculate_analyticalFCPCC(empSC, c_best, sigma_noise);
corr_max = corr(empPCC(tril(ones(number_of_nodes),-1)>0), aPCC(tril(ones(number_of_nodes),-1)>0)); % compute higher correlation removing the diagonal

end