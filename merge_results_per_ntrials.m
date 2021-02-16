function merge_results_per_ntrials(cfg_simulations, root_dir, nexperiments, signific_thresh)
%MERGE_RESULTS_PER_NTRIALS
%   Merges results from the computations on varying number of trials in
%   total.
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

%%% Compute p-values for effect present

pvalues_MI_effectpresent = nan(nexperiments, length(cfg_simulations.levels_nTrials));
pvalues_POS_effectpresent = nan(nexperiments, length(cfg_simulations.levels_nTrials));
pvalues_watson_effectpresent = nan(nexperiments, length(cfg_simulations.levels_nTrials));
pvalues_logregress_effectpresent = nan(nexperiments, length(cfg_simulations.levels_nTrials));

dispstat('init');
for iexp=1:nexperiments
    
    dispstat(sprintf('\nExperiment %d%',iexp));

    if isfile(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']))
        
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectpresent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectpresent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectpresent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectpresent');
        
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectpresent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectpresent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectpresent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectpresent');

        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr_effectpresent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_surr_distr_effectpresent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr_effectpresent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr_effectpresent');

        for iTrials=1:length(cfg_simulations.levels_nTrials)
            
            [~, pvalues_MI_effectpresent(iexp, iTrials)] = ttest(MIs_empirical_effectpresent(:, iTrials), MIs_chance_levels_effectpresent(:, iTrials), 'tail', 'right');
            [~, pvalues_POS_effectpresent(iexp, iTrials)] = ttest(POS_empirical_effectpresent(:, iTrials), POS_chance_levels_effectpresent(:, iTrials), 'tail', 'right');
            [~, pvalues_watson_effectpresent(iexp, iTrials)] = ttest(U2watson_empirical_effectpresent(:, iTrials), U2watson_chance_levels_effectpresent(:, iTrials), 'tail', 'right');
            [~, pvalues_logregress_effectpresent(iexp, iTrials)] = ttest(rms_logregress_empirical_effectpresent(:, iTrials), rms_logregress_chance_levels_effectpresent(:, iTrials), 'tail', 'right');
            
        end
    end
end

%%% Compute p-values for effect absent

pvalues_MI_effectabsent = nan(nexperiments, length(cfg_simulations.levels_nTrials));
pvalues_POS_effectabsent = nan(nexperiments, length(cfg_simulations.levels_nTrials));
pvalues_watson_effectabsent = nan(nexperiments, length(cfg_simulations.levels_nTrials));
pvalues_logregress_effectabsent = nan(nexperiments, length(cfg_simulations.levels_nTrials));

dispstat('init');
for iexp=1:nexperiments
    
    dispstat(sprintf('\nIteration %d%',iexp));
    
    if isfile(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']))
        
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectabsent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectabsent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectabsent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectabsent');
        
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectabsent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectabsent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectabsent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectabsent');

        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr_effectabsent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_surr_distr_effectabsent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr_effectabsent');
        load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr_effectabsent');
        
        for iTrials=1:length(cfg_simulations.levels_nTrials)
            
            [~, pvalues_MI_effectabsent(iexp, iTrials)] = ttest(MIs_empirical_effectabsent(:, iTrials), MIs_chance_levels_effectabsent(:, iTrials), 'tail', 'right');
            [~, pvalues_POS_effectabsent(iexp, iTrials)] = ttest(POS_empirical_effectabsent(:, iTrials), POS_chance_levels_effectabsent(:, iTrials), 'tail', 'right');
            [~, pvalues_watson_effectabsent(iexp, iTrials)] = ttest(U2watson_empirical_effectabsent(:, iTrials), U2watson_chance_levels_effectabsent(:, iTrials), 'tail', 'right');
            [~, pvalues_logregress_effectabsent(iexp, iTrials)] = ttest(rms_logregress_empirical_effectabsent(:, iTrials), rms_logregress_chance_levels_effectabsent(:, iTrials), 'tail', 'right');
            
        end
    end
end

%%% Compute False Positive rate

FPrate_ntrials_MI = nan(1, length(cfg_simulations.levels_nTrials));
FPrate_ntrials_POS = nan(1, length(cfg_simulations.levels_nTrials));
FPrate_ntrials_watson = nan(1, length(cfg_simulations.levels_nTrials));
FPrate_ntrials_logregress = nan(1, length(cfg_simulations.levels_nTrials));

for iTrials=1:length(cfg_simulations.levels_nTrials)

    FPrate_ntrials_MI(iTrials) = length(find(pvalues_MI_effectabsent(:, iTrials)<signific_thresh))/length(find(~isnan(pvalues_MI_effectabsent(:, iTrials))));
    FPrate_ntrials_POS(iTrials) = length(find(pvalues_POS_effectabsent(:, iTrials)<signific_thresh))/length(find(~isnan(pvalues_POS_effectabsent(:, iTrials))));
    FPrate_ntrials_watson(iTrials) = length(find(pvalues_watson_effectabsent(:, iTrials)<signific_thresh))/length(find(~isnan(pvalues_watson_effectabsent(:, iTrials))));
    FPrate_ntrials_logregress(iTrials) = length(find(pvalues_logregress_effectabsent(:, iTrials)<signific_thresh))/length(find(~isnan(pvalues_logregress_effectabsent(:, iTrials))));

end

if ~isdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_FPrate_ntrials_MI_phase_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'FPrate_ntrials_MI');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_FPrate_ntrials_POS_phase_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'FPrate_ntrials_POS');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_FPrate_ntrials_U2watson_phase_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'FPrate_ntrials_watson');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_FPrate_ntrials_rms_logregress_phase_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'FPrate_ntrials_logregress');

%%% Compute sensitivity

sensitivity_ntrials_MI = nan(1, length(cfg_simulations.levels_nTrials));
sensitivity_ntrials_POS = nan(1, length(cfg_simulations.levels_nTrials));
sensitivity_ntrials_watson = nan(1, length(cfg_simulations.levels_nTrials));
sensitivity_ntrials_logregress = nan(1, length(cfg_simulations.levels_nTrials));

for iTrials=1:length(cfg_simulations.levels_nTrials)

    sensitivity_ntrials_MI(iTrials) = length(find(pvalues_MI_effectpresent(:, iTrials)<signific_thresh))/length(find(~isnan(pvalues_MI_effectpresent(:, iTrials))));
    sensitivity_ntrials_POS(iTrials) = length(find(pvalues_POS_effectpresent(:, iTrials)<signific_thresh))/length(find(~isnan(pvalues_POS_effectpresent(:, iTrials))));
    sensitivity_ntrials_watson(iTrials) = length(find(pvalues_watson_effectpresent(:, iTrials)<signific_thresh))/length(find(~isnan(pvalues_watson_effectpresent(:, iTrials))));
    sensitivity_ntrials_logregress(iTrials) = length(find(pvalues_logregress_effectpresent(:, iTrials)<signific_thresh))/length(find(~isnan(pvalues_logregress_effectpresent(:, iTrials))));

end

if ~isdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_sensitivity_ntrials_MI_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'sensitivity_ntrials_MI');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_sensitivity_ntrials_POS_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'sensitivity_ntrials_POS');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_sensitivity_ntrials_U2watson_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'sensitivity_ntrials_watson');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_sensitivity_ntrials_rms_logregress_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'sensitivity_ntrials_logregress');

end
