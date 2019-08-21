function run_ctrl_heterogeneous_sigma(this_type)
% This function makes use of a modified version of Saggio's method -- which 
% derives the optimal value of c (global coupling strength), based on the
% maximum correlation  between eFc and aFc. In this work, ** we use different 
% values sigma for the noise amplitudes for ALL, HUB, or PERIPHERY nodes. **
% 
% The input argument this_type is an integer value for:
%     1: variable noise in all nodes
%     2: variable noise in hub nodes
%     3: variable noise in periphery nodes
% Acronyms: 
%          eSC: empirical Structural Connectivity
%          aSC: analytical Structural Connectivity
%          eFC: empirial Functional Connectivity - Pearson Correlation
%               Coefficients
%          aFC: analytical Functional Connectivity - Pearson Correlation
%               Coefficients
% NOTE: approx runtime 12 min on a 14 cores CPU using parfors
% =============================== CONTROLS ===============================%
% 1 - Load data for controls
seeds = [42, 51, 48, 91, 61, 62, 86, 81, 58, 19, 24, 89, 33, 49, 17];
for that_seed = 1:length(seeds)
    % 1 - Load data for adhd -- assumes the current working directory is
    % 'ADHDSCFC/Scripts/Model/'
    path_to_input_files = '../../Data/';
    
    load([path_to_input_files 'Schaeffer214-SC-FC/FC/CTRL-FC-PCC-118.mat'])
    load([path_to_input_files 'Schaeffer214-SC-FC/SC/CTRLSC118.mat'])
    load([path_to_input_files 'Schaeffer214-NodeEdgeClassification/ADHDSCConns.mat'])
    load([path_to_input_files 'Schaeffer214-NodeEdgeClassification/ADHDSCHubs.mat'])

    % Select node type
    node_types = {'all', 'hubs', 'periphery'};

    % Structural 
    SC    = CTRLSC; 
    clear CTRLSC

    % Functional Pearson Correlation Coefficients
    FCPCC = CTRL_FC_PCC; 
    clear CTRL_FC_PCC

    % Edge type
    HBL = hubMat.CTRL;

    % Node type
    HBN = hublist.CTRL;

        if strcmp('periphery', node_types{this_type})
            HBN = ~HBN;
        elseif strcmp('all', node_types{this_type})
            HBN = [];
        end

    % 2 - Load analytic SC -- asummed that it has been precalculateed
    load([path_to_input_files 'Schaeffer214-SC-FC/SC/CTRLaSC118.mat'])


    % 3 - Calculate analytic FC introducing variability in ALL nodes

    % Values of standard deviation of the distribution we use to introduce
    % variability in sigma of the noise;
    dist_std_values = linspace(0, 0.5, 32);

    % define output structures
    r_ctrl_esc_afc_hubs = struct();
    r_ctrl_asc_afc_hubs = struct();

    % preallocate memory for the matrices
    aCTRL_FC_HUBS = zeros([size(SC), length(dist_std_values)]);
    hubs_std_values = zeros([size(SC, 3), length(dist_std_values)]);

    % seed for getting random number for sigma
    this_seed = seeds(that_seed);

    % mean of the distriution of sigma_i
    mu = 1;

    for kk=1:length(dist_std_values)
        sigma_noise    = get_sigma_noise(size(SC, 1), dist_std_values(kk), this_seed, mu);

        [aCTRL_FC_HUBS(:, :, :, kk), c_best, c_critic, corr_max, c_range, correl, hubs_std_values(:, kk)] = perform_pse_subjects(SC, FCPCC, sigma_noise, HBN);

        %Calculate r_SC_FC for eSC-aFC with different levels of variability
        %Calculate r_SC_FC for aSC-aFC with different levels of variability
        for subj_idx=1:size(SC, 3)
            % eSC-aFC
            r_ctrl_esc_afc_hubs.all_edges(subj_idx, kk)       = calculate_corr_sc_fc(SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx, kk));
            r_ctrl_esc_afc_hubs.hub_edges(subj_idx, kk)       = calculate_corr_sc_fc(SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx, kk), HBL(:, :, :, subj_idx), 'hub');
            r_ctrl_esc_afc_hubs.feed_edges(subj_idx, kk)      = calculate_corr_sc_fc(SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx, kk), HBL(:, :, :, subj_idx), 'feeders');
            r_ctrl_esc_afc_hubs.periphery_edges(subj_idx, kk) = calculate_corr_sc_fc(SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx, kk), HBL(:, :, :, subj_idx), 'periphery');

            % aSC-aFC
            r_ctrl_asc_afc_hubs.all_edges(subj_idx, kk)       = calculate_corr_sc_fc(aCTRL_SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx, kk));
            r_ctrl_asc_afc_hubs.hub_edges(subj_idx, kk)       = calculate_corr_sc_fc(aCTRL_SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx, kk), HBL(:, :, :, subj_idx), 'hub');
            r_ctrl_asc_afc_hubs.feed_edges(subj_idx, kk)      = calculate_corr_sc_fc(aCTRL_SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx, kk), HBL(:, :, :, subj_idx), 'feeders');
            r_ctrl_asc_afc_hubs.periphery_edges(subj_idx, kk) = calculate_corr_sc_fc(aCTRL_SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx, kk), HBL(:, :, :, subj_idx), 'periphery');

        end
    end

    path_to_output_files = '../../Results/Schaeffer214-Model/';

    if strcmp('periphery', node_types{this_type})
        % Swap variable names and clean up after ourselves
        r_ctrl_asc_afc_periphery = r_ctrl_asc_afc_hubs; clear r_ctrl_asc_afc_hubs
        r_ctrl_esc_afc_periphery = r_ctrl_esc_afc_hubs; clear r_ctrl_esc_afc_hubs
        periphery_std_values = hubs_std_values; clear hubs_std_values

        % 4 - Save results
        filename = [path_to_output_files, 'CTRL_variable_noise_periphery_nodes_seed_', num2str(this_seed)];
        save(filename, 'dist_std_values', ...
                       'periphery_std_values', ...
                       'r_ctrl_asc_afc_periphery', ...
                       'r_ctrl_esc_afc_periphery');
                   
    elseif strcmp('hubs', node_types{this_type})
        % 4 - Save results
        filename = [path_to_output_files, 'CTRL_variable_noise_hub_nodes_seed_', num2str(this_seed)];
        save(filename, 'dist_std_values', ...
                       'hubs_std_values', ...
                       'r_ctrl_asc_afc_hubs', ...
                       'r_ctrl_esc_afc_hubs');
    else           
        % Swap variable names and clean up after ourselves
        r_ctrl_asc_afc = r_ctrl_asc_afc_hubs; clear r_ctrl_asc_afc_hubs
        r_ctrl_esc_afc = r_ctrl_esc_afc_hubs; clear r_ctrl_esc_afc_hubs
        all_std_values = hubs_std_values; clear hubs_std_values

        % 4 - Save results
        filename = [path_to_output_files,'CTRL_variable_noise_all_nodes_seed_', num2str(this_seed)];
        save(filename, 'dist_std_values', ...
                       'all_std_values', ...
                       'r_ctrl_asc_afc', ...
                       'r_ctrl_esc_afc');
    end
end
end % function run_ctrl_heterogeneous_sigma()
