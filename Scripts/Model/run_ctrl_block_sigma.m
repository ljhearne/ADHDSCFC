function run_ctrl_block_sigma()
% This function uses two values of sigma_i: one for the hubs and one for
% the periphery
% 
% The input argument is:
% Acronyms: 
%          eSC: empirical Structural Connectivity
%          aSC: analytical Structural Connectivity
%          eFC: empirial Functional Connectivity - Pearson Correlation
%               Coefficients
%          aFC: analytical Functional Connectivity - Pearson Correlation
%               Coefficients
% =============================== CONTROLS ===============================%
% 1 - Load data for controls
seeds = [42];  %51, 48, 91, 61, 62, 86, 81, 58, 19, 24, 89, 33, 49, 17];

for that_seed = 1:length(seeds)
    % 1 - Load data for adhd -- assumes the current working directory is
    % 'ADHDSCFC/Scripts/Model/'
    path_to_input_files = '../../Data/';
    
    load([path_to_input_files 'Schaeffer214-SC-FC/FC/CTRL-FC-PCC-118.mat'])
    load([path_to_input_files 'Schaeffer214-SC-FC/SC/CTRLSC118.mat'])
    load([path_to_input_files 'Schaeffer214-NodeEdgeClassification/ADHDSCConns.mat'])
    load([path_to_input_files 'Schaeffer214-NodeEdgeClassification/ADHDSCHubs.mat'])

    
    % Structural 
    SC    = CTRLSC; 
    clear CTRLSC

    % Functional Pearson Correlation Coefficients
    FCPCC = CTRL_FC_PCC; 
    clear CTRL_FC_PCC

    % Edge type
    HBL = hubMat.CTRL;

    % Get hub list 
    HBN = hublist.CTRL;

    % 2 - Load analytic SC -- asummed that it has been precalculateed
    load([path_to_input_files 'Schaeffer214-SC-FC/SC/CTRLaSC118.mat'])


    % 3 - Calculate analytic FC 

    % The standard deviation of the distribution of sigma_i is zero in this case 
    dist_std_value = 0; 
    
    % Generate meshgrid with values of sigma_h and sigma_p
    min_sigma = 2^-3;
    max_sigma = 2;
    
    % SH - sigma hubs / SP -- sigma periphery
    [SH, SP] = meshgrid(min_sigma:min_sigma:max_sigma, min_sigma:min_sigma:max_sigma);

    % define output structures
    r_ctrl_esc_afc_hubs = struct();
    r_ctrl_asc_afc_hubs = struct();

    % preallocate memory for the matrices
    aCTRL_FC_HUBS = zeros([size(SC), numel(SH)]);
    hubs_std_values = zeros([size(SC, 3), numel(SH)]);

    % seed for getting random number for sigma
    this_seed = seeds(that_seed);



    for kk=1:numel(SH)
        
        % All nodes get this value but in this function, only the hubs will get this value 
        mu_a = SH(kk);
        % Mean of the distribution of sigma_ for the periphery 
        mu_b = SP(kk);

        sigma_noise    = get_sigma_noise(size(SC, 1), dist_std_value, this_seed, mu_a);


        [aCTRL_FC_HUBS(:, :, :, kk), c_best, c_critic, corr_max, c_range, correl, hubs_std_values(:, kk)] = perform_pse_subjects(SC, FCPCC, sigma_noise, HBN, mu_b);

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

    % Swap variable names and clean up after ourselves
    r_ctrl_asc_afc = r_ctrl_asc_afc_hubs; clear r_ctrl_asc_afc_hubs
    r_ctrl_esc_afc = r_ctrl_esc_afc_hubs; clear r_ctrl_esc_afc_hubs
    all_std_values = hubs_std_values; clear hubs_std_values

    % 4 - Save results
    filename = [path_to_output_files,'CTRL_block_noise_seed_', num2str(this_seed)];
    save(filename, 'SH', 'SP', ...
                   'all_std_values', ...
                   'r_ctrl_asc_afc', ...
                   'r_ctrl_esc_afc');
    end

end % function run_ctrl_block_sigma()
