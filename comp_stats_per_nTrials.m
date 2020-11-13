function comp_stats_per_nTrials(cfg_simulations, iexp, phase_all_subjects, root_dir, ISI_range)
%COMP_STATS_PER_NTRIALS
%   Runs simulations to compare the sensitivity for different circular 
%   statistical tests (Phase Opposition Sum (POS), Circular Logistic 
%   Regression, Watson test, Modulation Index (MI) and Rayleigh test) as a
%   function of the number of trials in total.
%   Creates time series of hits and misses for each number of trials and 
%   computes empirical statistics and null distributions subject by 
%   subject. This is done with an injected effect and without injected
%   effect (hits and misses assigned purely randomly independently of
%   phase), to estimate sensitivity and False Positive rate respectively.
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

fprintf(['\n##### Experiment ' num2str(iexp) '\n']);

% Note number of subjects in total
nsubjects = length(phase_all_subjects);

% Initialize arrays and cells for empirical and surrogate phase-outcome
% statistics
% Effect present (used for estimating sensitivity)
MIs_empirical_effectpresent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
POS_empirical_effectpresent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
z_rayleigh_empirical_effectpresent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
U2watson_empirical_effectpresent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
rms_logregress_empirical_effectpresent = nan(nsubjects,length(cfg_simulations.levels_nTrials));

MIs_surr_distr_effectpresent = cell(nsubjects,length(cfg_simulations.levels_nTrials));
POS_surr_distr_effectpresent = cell(nsubjects,length(cfg_simulations.levels_nTrials));
z_rayleigh_surr_distr_effectpresent = cell(nsubjects,length(cfg_simulations.levels_nTrials));
U2watson_surr_distr_effectpresent = cell(nsubjects,length(cfg_simulations.levels_nTrials));
rms_logregress_surr_distr_effectpresent = cell(nsubjects,length(cfg_simulations.levels_nTrials));

MIs_chance_levels_effectpresent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
POS_chance_levels_effectpresent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
z_rayleigh_chance_levels_effectpresent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
U2watson_chance_levels_effectpresent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
rms_logregress_chance_levels_effectpresent = nan(nsubjects,length(cfg_simulations.levels_nTrials));

% Effect absent (used for estimating False Positive rate)
MIs_empirical_effectabsent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
POS_empirical_effectabsent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
z_rayleigh_empirical_effectabsent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
U2watson_empirical_effectabsent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
rms_logregress_empirical_effectabsent = nan(nsubjects,length(cfg_simulations.levels_nTrials));

MIs_surr_distr_effectabsent = cell(nsubjects,length(cfg_simulations.levels_nTrials));
POS_surr_distr_effectabsent = cell(nsubjects,length(cfg_simulations.levels_nTrials));
z_rayleigh_surr_distr_effectabsent = cell(nsubjects,length(cfg_simulations.levels_nTrials));
U2watson_surr_distr_effectabsent = cell(nsubjects,length(cfg_simulations.levels_nTrials));
rms_logregress_surr_distr_effectabsent = cell(nsubjects,length(cfg_simulations.levels_nTrials));

MIs_chance_levels_effectabsent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
POS_chance_levels_effectabsent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
z_rayleigh_chance_levels_effectabsent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
U2watson_chance_levels_effectabsent = nan(nsubjects,length(cfg_simulations.levels_nTrials));
rms_logregress_chance_levels_effectabsent = nan(nsubjects,length(cfg_simulations.levels_nTrials));

% Select trials for each trial number, for both effect present and absent.
% Here make sure that for each increase in number of trials we also
% include the same set of trials for the previous trial number. For 
% example when using 100 trials, these 100 trials also contain the same
% 75 trials from the previous level.
for isubject=1:nsubjects
    
    trials_by_trial_number_effectpresent_hits = cell(1, length(cfg_simulations.levels_nTrials));
    trials_by_trial_number_effectpresent_misses = cell(1, length(cfg_simulations.levels_nTrials));

    trials_by_trial_number_effectabsent_hits = cell(1, length(cfg_simulations.levels_nTrials));
    trials_by_trial_number_effectabsent_misses = cell(1, length(cfg_simulations.levels_nTrials));
    
    % Create phases for hits and misses
    [phases_hits_all,phases_misses_all] = distribute_outcomes(cfg_simulations, ISI_range, phase_all_subjects{isubject}, cfg_simulations.poc_simulations_Ntrials, max(cfg_simulations.levels_nTrials));
    
    % Go in inverse order through the trials, starting with maximum number
    % of trials, selecting random trials, then selecting a smaller subset
    % of these trials etc.
    for iTrials=length(cfg_simulations.levels_nTrials):-1:1

        nTrials = cfg_simulations.levels_nTrials(iTrials);

        if iTrials==length(cfg_simulations.levels_nTrials)
            trials_by_trial_number_effectpresent_hits{iTrials} = datasample(phases_hits_all, round(nTrials/2), 'Replace',false);
            trials_by_trial_number_effectpresent_misses{iTrials} = datasample(phases_misses_all, round(nTrials/2), 'Replace',false);
        else
            trials_by_trial_number_effectpresent_hits{iTrials} = datasample(trials_by_trial_number_effectpresent_hits{iTrials+1}, round(nTrials/2), 'Replace',false);
            trials_by_trial_number_effectpresent_misses{iTrials} = datasample(trials_by_trial_number_effectpresent_misses{iTrials+1}, round(nTrials/2), 'Replace',false);
        end
        
        allphases = [trials_by_trial_number_effectpresent_hits{iTrials} trials_by_trial_number_effectpresent_misses{iTrials}];
        
        % For effect absent condition:
        % Abolish the effect by randomly reshuffling trial labels
        order = randperm(length(allphases));
        perm1 = order(1:length(trials_by_trial_number_effectpresent_hits{iTrials}));
        perm2 = order(length(trials_by_trial_number_effectpresent_hits{iTrials})+1:end);

        trials_by_trial_number_effectabsent_hits{iTrials} = allphases(perm1);
        trials_by_trial_number_effectabsent_misses{iTrials} = allphases(perm2);

    end

    if ~isdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Phases_groups' filesep cfg_simulations.oscillation filesep]))
        mkdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Phases_groups' filesep cfg_simulations.oscillation filesep]))
    end
    
    save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Phases_groups' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_phases_hits_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_Subject' num2str(isubject) '_experiment' num2str(iexp) '.mat']), 'trials_by_trial_number_effectpresent_hits');
    save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Phases_groups' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_phases_misses_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_Subject' num2str(isubject) '_experiment' num2str(iexp) '.mat']), 'trials_by_trial_number_effectpresent_misses');

    save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Phases_groups' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_phases_hits_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_Subject' num2str(isubject) '_experiment' num2str(iexp) '.mat']), 'trials_by_trial_number_effectabsent_hits');
    save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Phases_groups' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_phases_misses_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_Subject' num2str(isubject) '_experiment' num2str(iexp) '.mat']), 'trials_by_trial_number_effectabsent_misses');

    % Compute stats for each number of trials, for effect present and absent condition
    for effectabsent=0:1

        fprintf('Computing stats for each number of trials...\n');
        for iTrials=1:length(cfg_simulations.levels_nTrials)
            
            % take trials (as preselected before)
            if effectabsent
                phases_hits = trials_by_trial_number_effectabsent_hits{iTrials};
                phases_misses = trials_by_trial_number_effectabsent_misses{iTrials};
            else
                phases_hits = trials_by_trial_number_effectpresent_hits{iTrials};
                phases_misses = trials_by_trial_number_effectpresent_misses{iTrials};
            end

            % compute stats
            if effectabsent
                
                [MIs_empirical_effectabsent(isubject, iTrials), POS_empirical_effectabsent(isubject, iTrials), z_rayleigh_empirical_effectabsent(isubject, iTrials), U2watson_empirical_effectabsent(isubject, iTrials), rms_logregress_empirical_effectabsent(isubject, iTrials)] = calc_phase_statistics(phases_hits, phases_misses, cfg_simulations, 0);

                [MIs_surr_distr_effectabsent{isubject, iTrials}, POS_surr_distr_effectabsent{isubject, iTrials}, z_rayleigh_surr_distr_effectabsent{isubject, iTrials}, U2watson_surr_distr_effectabsent{isubject, iTrials}, rms_logregress_surr_distr_effectabsent{isubject, iTrials}] = calc_chance_level_statistics(phases_hits, phases_misses, cfg_simulations, 0);

                MIs_chance_levels_effectabsent(isubject, iTrials) = nanmean(MIs_surr_distr_effectabsent{isubject, iTrials});            % for the moment ignore nan due to unsampled phase bin, but should consider threshold on percentage of nan
                POS_chance_levels_effectabsent(isubject, iTrials) = mean(POS_surr_distr_effectabsent{isubject, iTrials});
                z_rayleigh_chance_levels_effectabsent(isubject, iTrials) = mean(z_rayleigh_surr_distr_effectabsent{isubject, iTrials});
                U2watson_chance_levels_effectabsent(isubject, iTrials) = mean(U2watson_surr_distr_effectabsent{isubject, iTrials});
                rms_logregress_chance_levels_effectabsent(isubject, iTrials) = mean(rms_logregress_surr_distr_effectabsent{isubject, iTrials});

            else

                [MIs_empirical_effectpresent(isubject, iTrials), POS_empirical_effectpresent(isubject, iTrials), z_rayleigh_empirical_effectpresent(isubject, iTrials), U2watson_empirical_effectpresent(isubject, iTrials), rms_logregress_empirical_effectpresent(isubject, iTrials)] = calc_phase_statistics(phases_hits, phases_misses, cfg_simulations, 0);

                [MIs_surr_distr_effectpresent{isubject, iTrials}, POS_surr_distr_effectpresent{isubject, iTrials}, z_rayleigh_surr_distr_effectpresent{isubject, iTrials}, U2watson_surr_distr_effectpresent{isubject, iTrials}, rms_logregress_surr_distr_effectpresent{isubject, iTrials}] = calc_chance_level_statistics(phases_hits,phases_misses,cfg_simulations,0);

                MIs_chance_levels_effectpresent(isubject, iTrials) = nanmean(MIs_surr_distr_effectpresent{isubject, iTrials});            % for the moment ignore nan due to unsampled phase bin, but should consider threshold on percentage of nan
                POS_chance_levels_effectpresent(isubject, iTrials) = mean(POS_surr_distr_effectpresent{isubject, iTrials});
                z_rayleigh_chance_levels_effectpresent(isubject, iTrials) = mean(z_rayleigh_surr_distr_effectpresent{isubject, iTrials});
                U2watson_chance_levels_effectpresent(isubject, iTrials) = mean(U2watson_surr_distr_effectpresent{isubject, iTrials});
                rms_logregress_chance_levels_effectpresent(isubject, iTrials) = mean(rms_logregress_surr_distr_effectpresent{isubject, iTrials});
            end
        % end of loop through number of trials
        end
    % end of loop through effect present vs. absent
    end
% end of subject loop
end

if ~isdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep]))
end
if ~isdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep]))
end
if ~isdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep]))
end

% effect absent
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_empirical_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'z_rayleigh_empirical_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectabsent');

save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_chance_levels_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'z_rayleigh_chance_levels_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectabsent');

save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_surr_distr_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_surr_distr_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'z_rayleigh_surr_distr_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr_effectabsent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_effectabsent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr_effectabsent');

% effect present
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_empirical_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_empirical_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_empirical_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'z_rayleigh_empirical_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_empirical_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical_effectpresent');

save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_chance_levels_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_chance_levels_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'z_rayleigh_chance_levels_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels_effectpresent');

save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'POS_surr_distr_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_surr_distr_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'z_rayleigh_surr_distr_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr_effectpresent');
save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_effectpresent_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr_effectpresent');

if iexp==1      % save config-file only once
    save(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep cfg_simulations.oscillation '_cfg_simulations_ntrials_phaseoutcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '_experiment' num2str(iexp) '.mat']), 'cfg_simulations');
end
end