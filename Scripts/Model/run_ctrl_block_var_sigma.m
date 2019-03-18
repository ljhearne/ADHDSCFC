function run_ctrl_block_var_sigma()
% This function uses two values of sigma_i: one for the hubs and one for
% the periphery. In addition it introduces variability in the hubs and the
% periphery from two different distributions.
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
    
    % Generate meshgrid with values of std[sigma_h] and std[sigma_p]
    min_std_sigma_step_h = 2^-5;
    min_std_sigma_step_p = 2^-5;
    
    min_std_sigma_h = 2^-5;
    min_std_sigma_p = 2^-5;
    %
    max_std_sigma_h = 0.5;
    max_std_sigma_p = 0.5;
    
    % VSH - var sigma hubs / VSP -- var sigma periphery
    [VSH, VSP] = meshgrid(min_std_sigma_h:min_std_sigma_step_h:max_std_sigma_h, ...
                          min_std_sigma_p:min_std_sigma_step_p:max_std_sigma_p);

    % define output structures
    r_ctrl_esc_afc_hubs = struct();
    r_ctrl_asc_afc_hubs = struct();

    % preallocate memory for the matrices
    %aCTRL_FC_HUBS = zeros([size(SC), numel(SH)]);
    hubs_var_values = zeros([size(SC, 3), numel(VSH)]);

    % seed for getting random number for sigma
    this_seed = seeds(that_seed);


    for kk=1:numel(VSH)
        
         %  CTRL Mean of the distribution of sigma_i for the hubs
         mu_h = 0.4;
         %  CTRL Mean of the distribution of sigma_i for the periphery 
         mu_p = 0.5;
        
        % The range is actually the standard deviation
        std_hubs = VSH(kk);
        std_peri = VSP(kk);
        rng(this_seed)

        [aCTRL_FC_HUBS, c_best, c_critic, corr_max, c_range, correl, hubs_var_values(:, kk)] = perform_pse_subjects_block(SC, FCPCC, HBN, mu_h, mu_p, std_hubs, std_peri);

        %Calculate r_SC_FC for eSC-aFC with different levels of variability
        %Calculate r_SC_FC for aSC-aFC with different levels of variability
        for subj_idx=1:size(SC, 3)
            % eSC-aFC
            r_ctrl_esc_afc_hubs.all_edges(subj_idx, kk)       = calculate_corr_sc_fc(SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx));
            r_ctrl_esc_afc_hubs.hub_edges(subj_idx, kk)       = calculate_corr_sc_fc(SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx), HBL(:, :, :, subj_idx), 'hub');
            r_ctrl_esc_afc_hubs.feed_edges(subj_idx, kk)      = calculate_corr_sc_fc(SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx), HBL(:, :, :, subj_idx), 'feeders');
            r_ctrl_esc_afc_hubs.periphery_edges(subj_idx, kk) = calculate_corr_sc_fc(SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx), HBL(:, :, :, subj_idx), 'periphery');

            % aSC-aFC
            r_ctrl_asc_afc_hubs.all_edges(subj_idx, kk)       = calculate_corr_sc_fc(aCTRL_SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx));
            r_ctrl_asc_afc_hubs.hub_edges(subj_idx, kk)       = calculate_corr_sc_fc(aCTRL_SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx), HBL(:, :, :, subj_idx), 'hub');
            r_ctrl_asc_afc_hubs.feed_edges(subj_idx, kk)      = calculate_corr_sc_fc(aCTRL_SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx), HBL(:, :, :, subj_idx), 'feeders');
            r_ctrl_asc_afc_hubs.periphery_edges(subj_idx, kk) = calculate_corr_sc_fc(aCTRL_SC(:, :, subj_idx), aCTRL_FC_HUBS(:, :, subj_idx), HBL(:, :, :, subj_idx), 'periphery');

        end
    end

    path_to_output_files = '../../Results/Schaeffer214-Model/';

    % Swap variable names and clean up after ourselves
    r_ctrl_asc_afc = r_ctrl_asc_afc_hubs; clear r_ctrl_asc_afc_hubs
    r_ctrl_esc_afc = r_ctrl_esc_afc_hubs; clear r_ctrl_esc_afc_hubs
    all_var_values = hubs_var_values; clear hubs_var_values

    % 4 - Save results
    filename = [path_to_output_files,'CTRL_block_var_noise_16x16_hub_periphery_20_percent'];
    save(filename, 'VSH', 'VSP', ...
                   'all_var_values', ...
                   'r_ctrl_asc_afc', ...
                   'r_ctrl_esc_afc');
end

end % function run_ctrl_var_block_sigma()
