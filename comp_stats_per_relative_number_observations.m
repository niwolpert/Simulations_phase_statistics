function comp_stats_per_relative_number_observations(cfg_simulations, iexp, phase_all_subjects, root_dir, ISI_range)
%COMP_STATS_PER_RELATIVE_NUMBER_OBSERVATIONS
%   Runs simulations to compare the sensitivity for different circular 
%   statistical tests (Phase Opposition Sum (POS), Circular Logistic 
%   Regression, Watson test and Modulation Index (MI)) as a function of the 
%   relative number of observations in the two conditions, with vs. without 
%   resampling (the Rayleigh test is not included here as it is only a one
%   sample test) and with effect present and effect absent (to compute
%   sensitivity and False Positive rate respectively). Here, group level
%   statistics are estimated by running t-tests on empirical vs. chance
%   level phase-outcome statistics.
%   Saves empirical phase-outcome statistics and surrogate distributions.
%   
%   INPUTS
%   - cfg_simulations:                 Configuration structure with
%                                      simulation parameters
%   - iexp:                            Number of virtual experiment
%   - phase_all_subjects:              1xN cell with the oscillatory phase
%                                      time series of all N subjects
%   - root_dir:                        Root directory where results will be
%                                      saved
%   - ISI_range:                       The lower and upper limit for
%                                      intervals between events
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

rng('shuffle');

nsubjects = length(phase_all_subjects);

%%% Initialize results matrices

% effect present

MIs_empirical_effectpresent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_empirical_effectpresent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_empirical_effectpresent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_empirical_effectpresent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));

MIs_chance_levels_effectpresent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_chance_levels_effectpresent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_chance_levels_effectpresent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_chance_levels_effectpresent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));

MIs_surr_distr_effectpresent_with_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_surr_distr_effectpresent_with_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_surr_distr_effectpresent_with_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_surr_distr_effectpresent_with_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));

MIs_empirical_effectpresent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_empirical_effectpresent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_empirical_effectpresent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_empirical_effectpresent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));

MIs_chance_levels_effectpresent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_chance_levels_effectpresent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_chance_levels_effectpresent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_chance_levels_effectpresent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));

MIs_surr_distr_effectpresent_without_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_surr_distr_effectpresent_without_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_surr_distr_effectpresent_without_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_surr_distr_effectpresent_without_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));

% effect absent

MIs_empirical_effectabsent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_empirical_effectabsent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_empirical_effectabsent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_empirical_effectabsent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));

MIs_chance_levels_effectabsent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_chance_levels_effectabsent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_chance_levels_effectabsent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_chance_levels_effectabsent_with_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));

MIs_surr_distr_effectabsent_with_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_surr_distr_effectabsent_with_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_surr_distr_effectabsent_with_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_surr_distr_effectabsent_with_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));

MIs_empirical_effectabsent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_empirical_effectabsent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_empirical_effectabsent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_empirical_effectabsent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));

MIs_chance_levels_effectabsent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_chance_levels_effectabsent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_chance_levels_effectabsent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_chance_levels_effectabsent_without_resampling = nan(nsubjects, length(cfg_simulations.proportion_group_sizes));

MIs_surr_distr_effectabsent_without_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
POS_surr_distr_effectabsent_without_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
U2watson_surr_distr_effectabsent_without_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));
rms_logregress_surr_distr_effectabsent_without_resampling = cell(nsubjects, length(cfg_simulations.proportion_group_sizes));

for effectabsent=0:1
    
    for isubject=1:nsubjects
        
        % Create phases for hits and misses. Here, we first create a set of 
        % 240 hits and misses, and later subsample sets of hits and misses
        % to vary the relative number of observations
        if effectabsent
            phase_outcome_coupling = 0;
        else
            phase_outcome_coupling = cfg_simulations.poc_relative_trial_number;
        end
        [phases_hits_all,phases_misses_all] = distribute_outcomes(cfg_simulations, ISI_range, phase_all_subjects{isubject}, phase_outcome_coupling, (cfg_simulations.proportion_group_sizes(end)*cfg_simulations.nTrials_total_relative_trial_number)*2);
        
        % Prepare phases for each balance in relative number of trials.
        % From the pre-created sets of hits and misses, subsample subsets
        % such as to systematically vary the relative number of
        % observations for hits and misses, while keeping the total number
        % of trials constant (as specified in
        % cfg_simulations.nTrials_total_relative_trial_number).
        % Here make sure that for each increase in number of trials we also
        % include the same set of trials for the previous trial number. For 
        % example when using 100 trials, these 100 trials also contain the same
        % 75 trials.
        trials_by_trial_ratio_hits = cell(1, length(cfg_simulations.proportion_group_sizes));
        trials_by_trial_ratio_misses = cell(1, length(cfg_simulations.proportion_group_sizes));
        for i=length(cfg_simulations.proportion_group_sizes):-1:1
            
            nTrials = cfg_simulations.proportion_group_sizes(i)*cfg_simulations.nTrials_total_relative_trial_number;
            
            if i==length(cfg_simulations.proportion_group_sizes)
                
                trials_by_trial_ratio_hits{i} = datasample(phases_hits_all, nTrials, 'Replace',false);
                trials_by_trial_ratio_misses{i} = datasample(phases_misses_all, nTrials, 'Replace',false);
                
            else
                trials_by_trial_ratio_hits{i} = datasample(trials_by_trial_ratio_hits{i+1}, nTrials, 'Replace',false);
                trials_by_trial_ratio_misses{i} = datasample(trials_by_trial_ratio_misses{i+1}, nTrials, 'Replace',false);
            end
        end
        trials_by_trial_ratio_misses = fliplr(trials_by_trial_ratio_misses);
        
        if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'Phases' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep cfg_simulations.oscillation filesep]))
            mkdir(strcat([root_dir 'Relative_trial_number' filesep 'Phases' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep cfg_simulations.oscillation filesep]))
        end
        save(strcat([root_dir 'Relative_trial_number' filesep 'Phases' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_trials_by_trial_ratio_hits_effectabsent_level' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'trials_by_trial_ratio_hits');
        save(strcat([root_dir 'Relative_trial_number' filesep 'Phases' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_trials_by_trial_ratio_misses_effectabsent_level' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'trials_by_trial_ratio_misses');
        
        % Do the same computations with and without resampling
        for with_resampling=0:1
            
            if with_resampling

                fprintf(['\n### Computing stats with resampling\n']);
                nresamples=cfg_simulations.nresamples_relative_trial_number;
            else
                fprintf(['\n### Computing stats without resampling\n']);
                nresamples=0;
            end
            
            for igroupratio=1:length(cfg_simulations.proportion_group_sizes)
                
                % take trials (as preselected before)
                phases_hits = trials_by_trial_ratio_hits{igroupratio};
                phases_misses = trials_by_trial_ratio_misses{igroupratio};

                if ~effectabsent && ~with_resampling
                    
                    [MIs_empirical_effectpresent_without_resampling(isubject, igroupratio), POS_empirical_effectpresent_without_resampling(isubject, igroupratio), ~, U2watson_empirical_effectpresent_without_resampling(isubject, igroupratio), rms_logregress_empirical_effectpresent_without_resampling(isubject, igroupratio)] = calc_phase_statistics(phases_hits, phases_misses, cfg_simulations, nresamples);
                    
                    [MIs_surr_distr_effectpresent_without_resampling{isubject, igroupratio}, POS_surr_distr_effectpresent_without_resampling{isubject, igroupratio}, ~, U2watson_surr_distr_effectpresent_without_resampling{isubject, igroupratio}, rms_logregress_surr_distr_effectpresent_without_resampling{isubject, igroupratio}] = calc_chance_level_statistics(phases_hits, phases_misses, cfg_simulations, nresamples);
                    
                    MIs_chance_levels_effectpresent_without_resampling(isubject, igroupratio) = nanmean(MIs_surr_distr_effectpresent_without_resampling{isubject, igroupratio});
                    POS_chance_levels_effectpresent_without_resampling(isubject, igroupratio) = mean(POS_surr_distr_effectpresent_without_resampling{isubject, igroupratio});
                    U2watson_chance_levels_effectpresent_without_resampling(isubject, igroupratio) = mean(U2watson_surr_distr_effectpresent_without_resampling{isubject, igroupratio});
                    rms_logregress_chance_levels_effectpresent_without_resampling(isubject, igroupratio) = mean(rms_logregress_surr_distr_effectpresent_without_resampling{isubject, igroupratio});

                elseif ~effectabsent && with_resampling 
                    
                    [MIs_empirical_effectpresent_with_resampling(isubject, igroupratio), POS_empirical_effectpresent_with_resampling(isubject, igroupratio), ~, U2watson_empirical_effectpresent_with_resampling(isubject, igroupratio), rms_logregress_empirical_effectpresent_with_resampling(isubject, igroupratio)] = calc_phase_statistics(phases_hits, phases_misses, cfg_simulations, nresamples);
                    
                    [MIs_surr_distr_effectpresent_with_resampling{isubject, igroupratio}, POS_surr_distr_effectpresent_with_resampling{isubject, igroupratio}, ~, U2watson_surr_distr_effectpresent_with_resampling{isubject, igroupratio}, rms_logregress_surr_distr_effectpresent_with_resampling{isubject, igroupratio}] = calc_chance_level_statistics(phases_hits, phases_misses, cfg_simulations, nresamples);
                    
                    MIs_chance_levels_effectpresent_with_resampling(isubject, igroupratio) = nanmean(MIs_surr_distr_effectpresent_with_resampling{isubject, igroupratio});
                    POS_chance_levels_effectpresent_with_resampling(isubject, igroupratio) = mean(POS_surr_distr_effectpresent_with_resampling{isubject, igroupratio});
                    U2watson_chance_levels_effectpresent_with_resampling(isubject, igroupratio) = mean(U2watson_surr_distr_effectpresent_with_resampling{isubject, igroupratio});
                    rms_logregress_chance_levels_effectpresent_with_resampling(isubject, igroupratio) = mean(rms_logregress_surr_distr_effectpresent_with_resampling{isubject, igroupratio});
                    
                elseif effectabsent && ~with_resampling
                    
                    [MIs_empirical_effectabsent_without_resampling(isubject, igroupratio), POS_empirical_effectabsent_without_resampling(isubject, igroupratio), ~, U2watson_empirical_effectabsent_without_resampling(isubject, igroupratio), rms_logregress_empirical_effectabsent_without_resampling(isubject, igroupratio)] = calc_phase_statistics(phases_hits, phases_misses, cfg_simulations, nresamples);
                    
                    [MIs_surr_distr_effectabsent_without_resampling{isubject, igroupratio}, POS_surr_distr_effectabsent_without_resampling{isubject, igroupratio}, ~, U2watson_surr_distr_effectabsent_without_resampling{isubject, igroupratio}, rms_logregress_surr_distr_effectabsent_without_resampling{isubject, igroupratio}] = calc_chance_level_statistics(phases_hits, phases_misses, cfg_simulations, nresamples);
                    
                    MIs_chance_levels_effectabsent_without_resampling(isubject, igroupratio) = nanmean(MIs_surr_distr_effectabsent_without_resampling{isubject, igroupratio});
                    POS_chance_levels_effectabsent_without_resampling(isubject, igroupratio) = mean(POS_surr_distr_effectabsent_without_resampling{isubject, igroupratio});
                    U2watson_chance_levels_effectabsent_without_resampling(isubject, igroupratio) = mean(U2watson_surr_distr_effectabsent_without_resampling{isubject, igroupratio});
                    rms_logregress_chance_levels_effectabsent_without_resampling(isubject, igroupratio) = mean(rms_logregress_surr_distr_effectabsent_without_resampling{isubject, igroupratio});
                    
                elseif effectabsent && with_resampling
                    
                    [MIs_empirical_effectabsent_with_resampling(isubject, igroupratio), POS_empirical_effectabsent_with_resampling(isubject, igroupratio), ~, U2watson_empirical_effectabsent_with_resampling(isubject, igroupratio), rms_logregress_empirical_effectabsent_with_resampling(isubject, igroupratio)] = calc_phase_statistics(phases_hits, phases_misses, cfg_simulations, nresamples);
                    
                    [MIs_surr_distr_effectabsent_with_resampling{isubject, igroupratio}, POS_surr_distr_effectabsent_with_resampling{isubject, igroupratio}, ~, U2watson_surr_distr_effectabsent_with_resampling{isubject, igroupratio}, rms_logregress_surr_distr_effectabsent_with_resampling{isubject, igroupratio}] = calc_chance_level_statistics(phases_hits, phases_misses, cfg_simulations, nresamples);
                    
                    MIs_chance_levels_effectabsent_with_resampling(isubject, igroupratio) = nanmean(MIs_surr_distr_effectabsent_with_resampling{isubject, igroupratio});
                    POS_chance_levels_effectabsent_with_resampling(isubject, igroupratio) = mean(POS_surr_distr_effectabsent_with_resampling{isubject, igroupratio});
                    U2watson_chance_levels_effectabsent_with_resampling(isubject, igroupratio) = mean(U2watson_surr_distr_effectabsent_with_resampling{isubject, igroupratio});
                    rms_logregress_chance_levels_effectabsent_with_resampling(isubject, igroupratio) = mean(rms_logregress_surr_distr_effectabsent_with_resampling{isubject, igroupratio});
                    
                end
            end
        end
    end
end

%%% Effect present

% with resampling
if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep]))
end
if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep]))
end
if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep]))
end

save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectpresent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectpresent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectpresent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectpresent_with_resampling');

save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectpresent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectpresent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectpresent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectpresent_with_resampling');

save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr_effectpresent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_surr_distr_effectpresent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr_effectpresent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr_effectpresent_with_resampling');

% without resampling
if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep]))
end
if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep]))
end
if ~isdir(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep]))
end

save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectpresent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectpresent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectpresent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectpresent_without_resampling');

save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectpresent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectpresent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectpresent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectpresent_without_resampling');

save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr_effectpresent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_surr_distr_effectpresent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr_effectpresent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr_effectpresent_without_resampling');

%%% Effect absent

% with resampling
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectabsent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectabsent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectabsent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectabsent_with_resampling');

save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectabsent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectabsent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectabsent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectabsent_with_resampling');

save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr_effectabsent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_surr_distr_effectabsent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr_effectabsent_with_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_with_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr_effectabsent_with_resampling');

% without resampling
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectabsent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectabsent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectabsent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectabsent_without_resampling');

save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectabsent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectabsent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectabsent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectabsent_without_resampling');

save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr_effectabsent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'POS_surr_distr_effectabsent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_without_resampling_effectabsent__phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr_effectabsent_without_resampling');
save(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_without_resampling_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr_effectabsent_without_resampling');

if iexp==1      % save config-file only once
    save(strcat([root_dir 'Relative_trial_number' filesep cfg_simulations.oscillation '_cfg_simulations_relative_trial_number_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '_' num2str(nTrials_total_relative_trial_number) 'trialstotal.mat']), 'cfg_simulations');
end
end