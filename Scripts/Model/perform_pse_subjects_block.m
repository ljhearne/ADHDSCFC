function [aFCPCC, c_best, c_critic, corr_max, c_range, correl, var_sigmas] = perform_pse_subjects_block(SC, FCPCC, HBN, mu_a, mu_b, std_a, std_b)

% This is a dummy function to wrap the parfor and the parameter exploration
% to detect the optimal c-value across subjects when there is variability
% in the hubs and in the periphery

    number_of_subjects = size(SC, 3);
    number_of_nodes    = size(SC, 2);
    node_attributes = HBN;

    

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
    var_sigmas = zeros(number_of_subjects, 1);
    
    parfor subject_idx = 1:size(SC, 3)
        empSC   = SC(:, :, subject_idx);
        empPCC  = FCPCC(:, :, subject_idx);
        [Sigma(:, :, subject_idx), var_sigmas(subject_idx)] = introduce_block_heterogeneity(node_attributes(subject_idx, :), mu_a, mu_b, std_a, std_b);
        
        [aFCPCC(:, :, subject_idx), c_best(subject_idx), ...     % optimal value of coupling that maximizes similarity between empFC and anaFC
            c_critic(subject_idx), ...   % critical value of coupling derived from empSC
            corr_max(subject_idx), ...   % value of correlation between empFc and anaFC at c_best
            c_range(:, subject_idx), ... % range explored, there will be variability across subjects
            correl(:, subject_idx)] = perform_pse_coupling(empSC, empPCC, prec, Sigma(:, :, subject_idx));
        

    end
    clear Sigma
end