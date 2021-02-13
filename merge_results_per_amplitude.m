function merge_results_per_amplitude(cfg_simulations, root_dir, nexperiments, signific_thresh)
%MERGE_RESULTS_PER_AMPLITUDE
%   Merges results from the computations on varying amplitude of the
%   oscillation modulating outcomes
%   - Computes p-values on the group level by computing t-tests on 
%     empirical vs. chance level for each statistical test
%   - Computes sensitivity for each statistical test as the proportion of 
%     experiments with a phase-outcome coupling higher than zero yielding a 
%     significant p-value (<.05).
%   - Computes False Positive rate for each statistical test as the 
%     proportion of experiments yielding a significant p-value (<.05) in 
%     the absence of an injected effect (phase-outcome coupling = 0).
%   
%   INPUTS
%   - cfg_simulations:      Configuration structure with simulation 
%                           parameters
%   - root_dir:                        Root directory where results will be
%                                      saved
%   - nexperiments:         Number of virtual experiments
%   - signific_thresh:      Threshold for significance (default: 0.05)
%
% When using this function in any published study, please cite: Wolpert, 
% N., Tallon-Baudry, C. (2020). Evaluation of different statistical 
% procedures to estimate coupling between oscillatory phase and 
% behavioral response (in preparation)
%
% This function was written in Matlab version R2017b.
%
% Copyright (C) 2020, Laboratoire de Neurosciences Cognitives, Nicolai 
% Wolpert, Catherine Tallon-Baudry
% Email: nicolaiwolpert@gmail.com
% 
% DISCLAIMER:
% This code is provided without explicit or implicit guarantee, and without 
% any form of technical support. The code is not intended to be used for 
% clinical purposes. The functions are free to use and can be 
% redistributed, modified and adapted, under the terms of the CC BY-NC-SA
% version of creative commons license (see
% <https://creativecommons.org/licenses/>).

fprintf('\n### Merging files across iterations, computing fp-rate, sensitivity and d-prime...\n');

coupling_strengths_MI = nan(nexperiments, length(cfg_simulations.steps_amplitudes));
coupling_strengths_POS = nan(nexperiments, length(cfg_simulations.steps_amplitudes));
coupling_strengths_z_rayleigh = nan(nexperiments, length(cfg_simulations.steps_amplitudes));
coupling_strengths_U2watson = nan(nexperiments, length(cfg_simulations.steps_amplitudes));
coupling_strengths_rms_logregress = nan(nexperiments, length(cfg_simulations.steps_amplitudes));
coupling_strengths_KLD = nan(nexperiments, length(cfg_simulations.steps_amplitudes));

pvalues_coupling_strength_ttest_MI = nan(nexperiments, length(cfg_simulations.steps_amplitudes));
pvalues_coupling_strength_ttest_POS = nan(nexperiments, length(cfg_simulations.steps_amplitudes));
pvalues_coupling_strength_ttest_z_rayleigh = nan(nexperiments, length(cfg_simulations.steps_amplitudes));
pvalues_coupling_strength_ttest_U2watson = nan(nexperiments, length(cfg_simulations.steps_amplitudes));
pvalues_coupling_strength_ttest_rms_logregress = nan(nexperiments, length(cfg_simulations.steps_amplitudes));
pvalues_coupling_strength_ttest_KLD = nan(nexperiments, length(cfg_simulations.steps_amplitudes));

% get distribution of t-test t-values under pure noise
t_values_ttest_noise_MI = nan(1, nexperiments);
t_values_ttest_noise_POS = nan(1, nexperiments);
t_values_ttest_noise_rayleigh = nan(1, nexperiments);
t_values_ttest_noise_watson = nan(1, nexperiments);
t_values_ttest_noise_logregress = nan(1, nexperiments);
t_values_ttest_noise_KLD = nan(1, nexperiments);

dispstat('','init');
for iteration=1:nexperiments
    
    dispstat(sprintf('Iteration %d%',iteration));
    
    if isfile(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_MIs_empirical_iteration' num2str(iteration) '.mat']))
        
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_MIs_empirical_iteration' num2str(iteration) '.mat']), 'MIs_empirical');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_POS_empirical_iteration' num2str(iteration) '.mat']), 'POS_empirical');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_z_rayleigh_empirical_iteration' num2str(iteration) '.mat']), 'z_rayleigh_empirical');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_U2watson_empirical_iteration' num2str(iteration) '.mat']), 'U2watson_empirical');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_rms_logregress_empirical_iteration' num2str(iteration) '.mat']), 'rms_logregress_empirical');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_KLD_empirical_iteration' num2str(iteration) '.mat']), 'KLD_empirical');
        
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep 'amplitudes_MIs_chance_levels_iteration' num2str(iteration) '.mat']), 'MIs_chance_levels');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep 'amplitudes_POS_chance_levels_iteration' num2str(iteration) '.mat']), 'POS_chance_levels');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep 'amplitudes_z_rayleigh_chance_levels_iteration' num2str(iteration) '.mat']), 'z_rayleigh_chance_levels');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep 'amplitudes_U2watson_chance_levels_iteration' num2str(iteration) '.mat']), 'U2watson_chance_levels');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep 'amplitudes_rms_logregress_chance_levels_iteration' num2str(iteration) '.mat']), 'rms_logregress_chance_levels');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep 'amplitudes_KLD_chance_levels_iteration' num2str(iteration) '.mat']), 'KLD_chance_levels');
        
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep 'amplitudes_MIs_surr_distr_iteration' num2str(iteration) '.mat']), 'MIs_surr_distr');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep 'amplitudes_POS_surr_distr_iteration' num2str(iteration) '.mat']), 'POS_surr_distr');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep 'amplitudes_z_rayleigh_surr_distr_iteration' num2str(iteration) '.mat']), 'z_rayleigh_surr_distr');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep 'amplitudes_U2watson_surr_distr_iteration' num2str(iteration) '.mat']), 'U2watson_surr_distr');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep 'amplitudes_rms_logregress_surr_distr_iteration' num2str(iteration) '.mat']), 'rms_logregress_surr_distr');
        load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep 'amplitudes_KLD_surr_distr_iteration' num2str(iteration) '.mat']), 'KLD_surr_distr');
        
        for iSNR=1:length(cfg_simulations.steps_amplitudes)
            
            coupling_strengths_MI(iteration, iSNR) = mean(MIs_empirical(:, iSNR))-mean(MIs_chance_levels(:, iSNR));
            coupling_strengths_POS(iteration, iSNR) = mean(POS_empirical(:, iSNR))-mean(POS_chance_levels(:, iSNR));
            coupling_strengths_z_rayleigh(iteration, iSNR) = mean(z_rayleigh_empirical(:, iSNR))-mean(z_rayleigh_chance_levels(:, iSNR));
            coupling_strengths_U2watson(iteration, iSNR) = mean(U2watson_empirical(:, iSNR))-mean(U2watson_chance_levels(:, iSNR));
            coupling_strengths_rms_logregress(iteration, iSNR) = mean(rms_logregress_empirical(:, iSNR))-mean(rms_logregress_chance_levels(:, iSNR));
            coupling_strengths_KLD(iteration, iSNR) = mean(KLD_empirical(:, iSNR))-mean(KLD_chance_levels(:, iSNR));
            
            [~,pvalues_coupling_strength_ttest_MI(iteration, iSNR), ~, stat_tt_MI] = ttest(MIs_empirical(:, iSNR), MIs_chance_levels(:, iSNR),'Tail','right');
            [~,pvalues_coupling_strength_ttest_POS(iteration, iSNR), ~, stat_tt_POS] = ttest(POS_empirical(:, iSNR), POS_chance_levels(:, iSNR),'Tail','right');
            [~,pvalues_coupling_strength_ttest_z_rayleigh(iteration, iSNR), ~, stat_tt_rayleigh] = ttest(z_rayleigh_empirical(:, iSNR), z_rayleigh_chance_levels(:, iSNR),'Tail','right');
            [~,pvalues_coupling_strength_ttest_U2watson(iteration, iSNR), ~, stat_tt_watson] = ttest(U2watson_empirical(:, iSNR), U2watson_chance_levels(:, iSNR),'Tail','right');
            [~,pvalues_coupling_strength_ttest_rms_logregress(iteration, iSNR), ~, stat_tt_logregress] = ttest(rms_logregress_empirical(:, iSNR), rms_logregress_chance_levels(:, iSNR),'Tail','right');
            [~,pvalues_coupling_strength_ttest_KLD(iteration, iSNR), ~, stat_tt_KLD] = ttest(KLD_empirical(:, iSNR), KLD_chance_levels(:, iSNR),'Tail','right');
            
            % note t/z-values for pure noise
            if iSNR==length(cfg_simulations.steps_amplitudes)
                
                t_values_ttest_noise_MI(iteration) = stat_tt_MI.tstat;
                t_values_ttest_noise_POS(iteration) = stat_tt_POS.tstat;
                t_values_ttest_noise_rayleigh(iteration) = stat_tt_rayleigh.tstat;
                t_values_ttest_noise_watson(iteration) = stat_tt_watson.tstat;
                t_values_ttest_noise_logregress(iteration) = stat_tt_logregress.tstat;
                t_values_ttest_noise_KLD(iteration) = stat_tt_KLD.tstat;
                
            end
        end
    end
end

% compute sensitivity
sensitivity_coupling_strength_ttest_MI = nan(1, length(cfg_simulations.steps_amplitudes));
sensitivity_coupling_strength_ttest_POS = nan(1, length(cfg_simulations.steps_amplitudes));
sensitivity_coupling_strength_ttest_z_rayleigh = nan(1, length(cfg_simulations.steps_amplitudes));
sensitivity_coupling_strength_ttest_U2watson = nan(1, length(cfg_simulations.steps_amplitudes));
sensitivity_coupling_strength_ttest_rms_logregress = nan(1, length(cfg_simulations.steps_amplitudes));
sensitivity_coupling_strength_ttest_KLD = nan(1, length(cfg_simulations.steps_amplitudes));

for iSNR=1:length(cfg_simulations.steps_amplitudes)
    
    sensitivity_coupling_strength_ttest_MI(iSNR) = length(find(pvalues_coupling_strength_ttest_MI(:, iSNR)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_ttest_MI(:, iSNR))));
    sensitivity_coupling_strength_ttest_POS(iSNR) = length(find(pvalues_coupling_strength_ttest_POS(:, iSNR)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_ttest_POS(:, iSNR))));
    sensitivity_coupling_strength_ttest_z_rayleigh(iSNR) = length(find(pvalues_coupling_strength_ttest_z_rayleigh(:, iSNR)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_ttest_z_rayleigh(:, iSNR))));
    sensitivity_coupling_strength_ttest_U2watson(iSNR) = length(find(pvalues_coupling_strength_ttest_U2watson(:, iSNR)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_ttest_U2watson(:, iSNR))));
    sensitivity_coupling_strength_ttest_rms_logregress(iSNR) = length(find(pvalues_coupling_strength_ttest_rms_logregress(:, iSNR)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_ttest_rms_logregress(:, iSNR))));
    sensitivity_coupling_strength_ttest_KLD(iSNR) = length(find(pvalues_coupling_strength_ttest_KLD(:, iSNR)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_ttest_KLD(:, iSNR))));
    
end

%%% save

% p-values
if ~isdir(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep 'amplitudes_pvalues_coupling_strength_ttest_MI.mat']), 'pvalues_coupling_strength_ttest_MI');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep 'amplitudes_pvalues_coupling_strength_ttest_POS.mat']), 'pvalues_coupling_strength_ttest_POS');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep 'amplitudes_pvalues_coupling_strength_ttest_z_rayleigh.mat']), 'pvalues_coupling_strength_ttest_z_rayleigh');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep 'amplitudes_pvalues_coupling_strength_ttest_U2watson.mat']), 'pvalues_coupling_strength_ttest_U2watson');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep 'amplitudes_pvalues_coupling_strength_ttest_rms_logregress.mat']), 'pvalues_coupling_strength_ttest_rms_logregress');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep 'amplitudes_pvalues_coupling_strength_ttest_KLD.mat']), 'pvalues_coupling_strength_ttest_KLD');

% Sensitivity
if ~isdir(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep 'amplitudes_sensitivity_coupling_strength_ttest_MI.mat']), 'sensitivity_coupling_strength_ttest_MI');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep 'amplitudes_sensitivity_coupling_strength_ttest_POS.mat']), 'sensitivity_coupling_strength_ttest_POS');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep 'amplitudes_sensitivity_coupling_strength_ttest_z_rayleigh.mat']), 'sensitivity_coupling_strength_ttest_z_rayleigh');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep 'amplitudes_sensitivity_coupling_strength_ttest_U2watson.mat']), 'sensitivity_coupling_strength_ttest_U2watson');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep 'amplitudes_sensitivity_coupling_strength_ttest_rms_logregress.mat']), 'sensitivity_coupling_strength_ttest_rms_logregress');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep 'amplitudes_sensitivity_coupling_strength_ttest_KLD.mat']), 'sensitivity_coupling_strength_ttest_KLD');

end