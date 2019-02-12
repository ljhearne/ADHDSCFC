function [aFCPCC, c_best, c_critic, corr_max, c_range, correl, std_sigmas] = perform_pse_subjects(SC, FCPCC, sigma_noise, HBN, mu_b)

% This is a dummy function to wrap the parfor and the parameter exploration
% to detect the optimal c-value across subjects

    number_of_subjects = size(SC, 3);
    number_of_nodes    = size(SC, 2);
    
    if nargin < 2
        sigma_noise = 1; % uniform sigma for all nodes
    end
    
    if nargin < 4 || isempty(HBN)
        % Hack to make this function work for all cases
        node_attributes = ones(number_of_subjects, number_of_nodes);
        
    else
        node_attributes = HBN;
    end
    
    if nargin < 5
        % This value of mu is used for the nodes that do not have the
        % attribute specified in HBN
        mu_b = 1;
    end
    

    prec = 0.05;     % Percentage below c_critic at which the range to explore will be truncated

    % Preallocate output variables
    number_of_points = 1025; % NOTE: Hardoced value of points to explore between cmin and cmax.
    c_range = zeros(number_of_points, number_of_subjects);
    correl =  zeros(number_of_points, number_of_subjects);

    % Preallocate memory
    c_best   = zeros(number_of_subjects,1);
    c_critic = c_best;
    corr_max = c_best;
    aFCPCC   = zeros(size(FCPCC));

    % the code below should take approx 30s
    Sigma = zeros(size(SC));
    std_sigmas = zeros(number_of_subjects, 1);
    
    parfor subject_idx = 1:size(SC, 3)
        empSC   = SC(:, :, subject_idx);
        empPCC  = FCPCC(:, :, subject_idx);
        [Sigma(:, :, subject_idx), std_sigmas(subject_idx)] = introduce_node_based_variability(sigma_noise, node_attributes(subject_idx, :), mu_b);
        
        [aFCPCC(:, :, subject_idx), c_best(subject_idx), ...     % optimal value of coupling that maximizes similarity between empFC and anaFC
            c_critic(subject_idx), ...   % critical value of coupling derived from empSC
            corr_max(subject_idx), ...   % value of correlation between empFc and anaFC at c_best
            c_range(:, subject_idx), ... % range explored, there will be variability across subjects
            correl(:, subject_idx)] = perform_pse_coupling(empSC, empPCC, prec, Sigma(:, :, subject_idx));
        

    end
    clear Sigma
end