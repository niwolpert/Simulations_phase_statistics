function merge_results_per_relative_number_observations(cfg_simulations, root_dir, nexperiments, thresh_sign)
%MERGE_RESULTS_PER_RELATIVE_NUMBER_OBSERVATIONS
%   Merges results from the computations on varying relative number of
%   observations, with vs. without resampling.
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

fprintf('\nMerging files across subjects, computing sensitivity and FP-rate...\n');
pvalues_MI_effectpresent_withresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_POS_effectpresent_withresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_watson_effectpresent_withresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_logregress_effectpresent_withresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));

pvalues_MI_effectpresent_withoutresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_POS_effectpresent_withoutresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_watson_effectpresent_withoutresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_logregress_effectpresent_withoutresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));

pvalues_MI_effectabsent_withresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_POS_effectabsent_withresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_watson_effectabsent_withresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_logregress_effectabsent_withresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));

pvalues_MI_effectabsent_withoutresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_POS_effectabsent_withoutresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_watson_effectabsent_withoutresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));
pvalues_logregress_effectabsent_withoutresampling = nan(nexperiments, length(cfg_simulations.proportion_group_sizes));

dispstat('init');
for iexp=1:nexperiments
    
    dispstat(sprintf('\nExperiment %d%',iexp));
    
    if isfile(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_effectpresent_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']))
        
        % Effect present

        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectpresent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectpresent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectpresent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectpresent_with_resampling');
        
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectpresent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectpresent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectpresent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectpresent_with_resampling');
        
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectpresent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectpresent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectpresent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectpresent_without_resampling');

        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectpresent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectpresent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectpresent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectpresent_without_resampling');
        
        % Effect absent

        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectabsent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectabsent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectabsent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectabsent_with_resampling');

        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectabsent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectabsent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectabsent_with_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectabsent_with_resampling');
        
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectabsent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectabsent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectabsent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectabsent_without_resampling');

        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectabsent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectabsent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectabsent_without_resampling');
        load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectabsent_without_resampling');

        for igroupratio=1:length(cfg_simulations.proportion_group_sizes)
            
            [~, pvalues_MI_effectpresent_withresampling(iexp, igroupratio)] = ttest(MIs_empirical_effectpresent_with_resampling(:, igroupratio), MIs_chance_levels_effectpresent_with_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_POS_effectpresent_withresampling(iexp, igroupratio)] = ttest(POS_empirical_effectpresent_with_resampling(:, igroupratio), POS_chance_levels_effectpresent_with_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_watson_effectpresent_withresampling(iexp, igroupratio)] = ttest(U2watson_empirical_effectpresent_with_resampling(:, igroupratio), U2watson_chance_levels_effectpresent_with_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_logregress_effectpresent_withresampling(iexp, igroupratio)] = ttest(rms_logregress_empirical_effectpresent_with_resampling(:, igroupratio), rms_logregress_chance_levels_effectpresent_with_resampling(:, igroupratio), 'tail', 'right');
            
            [~, pvalues_MI_effectpresent_withoutresampling(iexp, igroupratio)] = ttest(MIs_empirical_effectpresent_without_resampling(:, igroupratio), MIs_chance_levels_effectpresent_without_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_POS_effectpresent_withoutresampling(iexp, igroupratio)] = ttest(POS_empirical_effectpresent_without_resampling(:, igroupratio), POS_chance_levels_effectpresent_without_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_watson_effectpresent_withoutresampling(iexp, igroupratio)] = ttest(U2watson_empirical_effectpresent_without_resampling(:, igroupratio), U2watson_chance_levels_effectpresent_without_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_logregress_effectpresent_withoutresampling(iexp, igroupratio)] = ttest(rms_logregress_empirical_effectpresent_without_resampling(:, igroupratio), rms_logregress_chance_levels_effectpresent_without_resampling(:, igroupratio), 'tail', 'right');
            
            [~, pvalues_MI_effectabsent_withresampling(iexp, igroupratio)] = ttest(MIs_empirical_effectabsent_with_resampling(:, igroupratio), MIs_chance_levels_effectabsent_with_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_POS_effectabsent_withresampling(iexp, igroupratio)] = ttest(POS_empirical_effectabsent_with_resampling(:, igroupratio), POS_chance_levels_effectabsent_with_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_watson_effectabsent_withresampling(iexp, igroupratio)] = ttest(U2watson_empirical_effectabsent_with_resampling(:, igroupratio), U2watson_chance_levels_effectabsent_with_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_logregress_effectabsent_withresampling(iexp, igroupratio)] = ttest(rms_logregress_empirical_effectabsent_with_resampling(:, igroupratio), rms_logregress_chance_levels_effectabsent_with_resampling(:, igroupratio), 'tail', 'right');
            
            [~, pvalues_MI_effectabsent_withoutresampling(iexp, igroupratio)] = ttest(MIs_empirical_effectabsent_without_resampling(:, igroupratio), MIs_chance_levels_effectabsent_without_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_POS_effectabsent_withoutresampling(iexp, igroupratio)] = ttest(POS_empirical_effectabsent_without_resampling(:, igroupratio), POS_chance_levels_effectabsent_without_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_watson_effectabsent_withoutresampling(iexp, igroupratio)] = ttest(U2watson_empirical_effectabsent_without_resampling(:, igroupratio), U2watson_chance_levels_effectabsent_without_resampling(:, igroupratio), 'tail', 'right');
            [~, pvalues_logregress_effectabsent_withoutresampling(iexp, igroupratio)] = ttest(rms_logregress_empirical_effectabsent_without_resampling(:, igroupratio), rms_logregress_chance_levels_effectabsent_without_resampling(:, igroupratio), 'tail', 'right');
            
        end
    end
end

FPrate_MI_withresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
FPrate_POS_withresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
FPrate_watson_withresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
FPrate_logregress_withresampling = nan(1, length(cfg_simulations.proportion_group_sizes));

sensitivity_MI_withresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
sensitivity_POS_withresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
sensitivity_watson_withresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
sensitivity_logregress_withresampling = nan(1, length(cfg_simulations.proportion_group_sizes));

FPrate_MI_withoutresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
FPrate_POS_withoutresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
FPrate_watson_withoutresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
FPrate_logregress_withoutresampling = nan(1, length(cfg_simulations.proportion_group_sizes));

sensitivity_MI_withoutresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
sensitivity_POS_withoutresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
sensitivity_watson_withoutresampling = nan(1, length(cfg_simulations.proportion_group_sizes));
sensitivity_logregress_withoutresampling = nan(1, length(cfg_simulations.proportion_group_sizes));

for igroupratio=1:length(cfg_simulations.proportion_group_sizes)
    
    sensitivity_MI_withresampling(igroupratio) = length(find(pvalues_MI_effectpresent_withresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_MI_effectpresent_withresampling(:, igroupratio))));
    sensitivity_POS_withresampling(igroupratio) = length(find(pvalues_POS_effectpresent_withresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_POS_effectpresent_withresampling(:, igroupratio))));
    sensitivity_watson_withresampling(igroupratio) = length(find(pvalues_watson_effectpresent_withresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_watson_effectpresent_withresampling(:, igroupratio))));
    sensitivity_logregress_withresampling(igroupratio) = length(find(pvalues_logregress_effectpresent_withresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_logregress_effectpresent_withresampling(:, igroupratio))));
    
    sensitivity_MI_withoutresampling(igroupratio) = length(find(pvalues_MI_effectpresent_withoutresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_MI_effectpresent_withoutresampling(:, igroupratio))));
    sensitivity_POS_withoutresampling(igroupratio) = length(find(pvalues_POS_effectpresent_withoutresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_POS_effectpresent_withoutresampling(:, igroupratio))));
    sensitivity_watson_withoutresampling(igroupratio) = length(find(pvalues_watson_effectpresent_withoutresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_watson_effectpresent_withoutresampling(:, igroupratio))));
    sensitivity_logregress_withoutresampling(igroupratio) = length(find(pvalues_logregress_effectpresent_withoutresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_logregress_effectpresent_withoutresampling(:, igroupratio))));
    
    FPrate_MI_withresampling(igroupratio) = length(find(pvalues_MI_effectabsent_withresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_MI_effectabsent_withresampling(:, igroupratio))));
    FPrate_POS_withresampling(igroupratio) = length(find(pvalues_POS_effectabsent_withresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_POS_effectabsent_withresampling(:, igroupratio))));
    FPrate_watson_withresampling(igroupratio) = length(find(pvalues_watson_effectabsent_withresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_watson_effectabsent_withresampling(:, igroupratio))));
    FPrate_logregress_withresampling(igroupratio) = length(find(pvalues_logregress_effectabsent_withresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_logregress_effectabsent_withresampling(:, igroupratio))));
    
    FPrate_MI_withoutresampling(igroupratio) = length(find(pvalues_MI_effectabsent_withoutresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_MI_effectabsent_withoutresampling(:, igroupratio))));
    FPrate_POS_withoutresampling(igroupratio) = length(find(pvalues_POS_effectabsent_withoutresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_POS_effectabsent_withoutresampling(:, igroupratio))));
    FPrate_watson_withoutresampling(igroupratio) = length(find(pvalues_watson_effectabsent_withoutresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_watson_effectabsent_withoutresampling(:, igroupratio))));
    FPrate_logregress_withoutresampling(igroupratio) = length(find(pvalues_logregress_effectabsent_withoutresampling(:, igroupratio)<thresh_sign))/length(find(~isnan(pvalues_logregress_effectabsent_withoutresampling(:, igroupratio))));
    
end

if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_MI_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_MI_withresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_POS_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_POS_withresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_watson_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_watson_withresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_logregress_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_logregress_withresampling');

if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_MI_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_MI_withoutresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_POS_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_POS_withoutresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_watson_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_watson_withoutresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_logregress_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_logregress_withoutresampling');

if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_MI_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_MI_withresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_POS_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_POS_withresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_watson_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_watson_withresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_logregress_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_logregress_withresampling');

if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_MI_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_MI_withoutresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_POS_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_POS_withoutresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_watson_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_watson_withoutresampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_logregress_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_logregress_withoutresampling');

end