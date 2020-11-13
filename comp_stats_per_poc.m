function comp_stats_per_poc(cfg_simulations, iexp, phase_all_subjects, root_dir, ISI_range)
%COMP_STATS_PER_POC
%   Runs simulations to compare the sensitivity for different circular 
%   statistical tests (Phase Opposition Sum (POS), Circular Logistic 
%   Regression, Watson test, Modulation Index (MI) and Rayleigh test).
%   Creates time series of hits and misses for each phase-outcome coupling
%   strength and computes empirical statistics and null distributions
%   subject by subject.
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

% Note number of subjects in total
nsubjects = length(phase_all_subjects);

% Initialize arrays and cells for empirical and surrogate phase-outcome
% statistics
MIs_empirical = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
POS_empirical = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
z_rayleigh_empirical = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
U2watson_empirical = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
rms_logregress_empirical = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));

MIs_surr_distr = cell(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
POS_surr_distr = cell(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
z_rayleigh_surr_distr = cell(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
U2watson_surr_distr = cell(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
rms_logregress_surr_distr = cell(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));

MIs_chance_levels = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
POS_chance_levels = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
z_rayleigh_chance_levels = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
U2watson_chance_levels = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
rms_logregress_chance_levels = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));

% p-values from direct output for rayleigh and watson
rayleigh_pvalues_direct = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
watson_pvalues_direct = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
logregress_pvalues_direct = nan(nsubjects, length(cfg_simulations.strengths_phase_outcome_coupling));
for isubject=1:nsubjects
    
    fprintf(['\n### Computing stats for Subject ' num2str(isubject) '\n']);
    
    % Create phases for hits and misses. Here, we start out with a 
    % phase-outcome coupling strength of 100 to systematically vary this
    % parameter in a next step
    [~,~,~,~,events_original] = distribute_outcomes(cfg_simulations, ISI_range, phase_all_subjects{isubject}, 100, cfg_simulations.nTrials_poc);
    
    ind_stim = find(events_original~=0);
    
    %%% Prepare sets of phases of hits and misses for each level of
    %%% phase-outcome coupling strength: 
    all_phases_hits_by_phaseoutcome_coupling = cell(1, length(cfg_simulations.strengths_phase_outcome_coupling));
    all_phases_misses_by_phaseoutcome_coupling = cell(1, length(cfg_simulations.strengths_phase_outcome_coupling));
    
    % To vary phase-outcome coupling strength, we systematically decrease
    % the percentage of 'noise trials', i.e. trials where outcome does not
    % depend on phase but is randomly assigned
    noise_trials = ind_stim;   % for first iteration, the noise trials are simply all the trials
    noise_trials_previous = noise_trials;
    for ipoc=1:length(cfg_simulations.strengths_phase_outcome_coupling)
        
        % compute the percentage and number of noise trials
        perc_noise_trials = 100-cfg_simulations.strengths_phase_outcome_coupling(ipoc);
        nNoiseTrials = round(cfg_simulations.nTrials_poc*(perc_noise_trials/100));
        
        % Imposing noise to reach the given phase-outcome coupling results
        % in differences in the number of trials for hits and misses. Here,
        % make sure that this difference in number of trials is not more
        % than 10%.
        percentage_diff_trial_number=Inf;
        while percentage_diff_trial_number>0.1
            
            events = events_original;
            
            % Select noise trials:
            % Take a subset of the noise trials for the previous noise level.
            % This way it is ensured that for every noise level iteration 
            % all the noise trials of the previous iteration with fewer noise
            % trials are included.
            noise_trials = sort(datasample(noise_trials_previous, nNoiseTrials, 'Replace',false));
            
            % Impose noise by switching trial outcome to hit/miss with 50% chance
            for iNoiseTrial = 1:nNoiseTrials

                events(noise_trials(iNoiseTrial)) = prob_select([0.5 0.5],[-1 1]);
                
            end
            
            ind_hits = find(events==1);
            ind_misses = find(events==-1);
            nhits=length(ind_hits);
            nmisses=length(ind_misses);
            
            diff_trial_number = abs(nhits - nmisses);
            percentage_diff_trial_number = diff_trial_number/max([nhits nmisses]);
        end
        
        all_phases_hits_by_phaseoutcome_coupling{ipoc} = phase_all_subjects{isubject}(ind_hits);
        all_phases_misses_by_phaseoutcome_coupling{ipoc} = phase_all_subjects{isubject}(ind_misses);
        noise_trials_previous = noise_trials;
        
    end
    
    % Save phases for hits and misses
    if ~isdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Phases_groups' filesep cfg_simulations.oscillation filesep]))
        mkdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Phases_groups' filesep cfg_simulations.oscillation filesep]))
    end
    save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Phases_groups' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_all_phases_hits_by_phaseoutcome_coupling_experiment' num2str(iexp) '_Subject' num2str(isubject) '.mat']), 'all_phases_hits_by_phaseoutcome_coupling');
    save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Phases_groups' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_all_phases_misses_by_phaseoutcome_coupling_experiment' num2str(iexp) '_Subject' num2str(isubject) '.mat']), 'all_phases_misses_by_phaseoutcome_coupling');
    
    % Compute stats for each phase-outcome coupling
    for ipoc=1:length(cfg_simulations.strengths_phase_outcome_coupling)
        
        phases_hits=all_phases_hits_by_phaseoutcome_coupling{ipoc};
        phases_misses=all_phases_misses_by_phaseoutcome_coupling{ipoc};
        
        [MIs_empirical(isubject, ipoc), POS_empirical(isubject, ipoc), z_rayleigh_empirical(isubject, ipoc), U2watson_empirical(isubject, ipoc), rms_logregress_empirical(isubject, ipoc), rayleigh_pvalues_direct(isubject, ipoc), watson_pvalues_direct(isubject, ipoc), logregress_pvalues_direct(isubject, ipoc)] = calc_phase_statistics(phases_hits, phases_misses, cfg_simulations, 0);
        
        [MIs_surr_distr{isubject, ipoc}, POS_surr_distr{isubject, ipoc}, z_rayleigh_surr_distr{isubject, ipoc}, U2watson_surr_distr{isubject, ipoc}, rms_logregress_surr_distr{isubject, ipoc}] = calc_chance_level_statistics(phases_hits,phases_misses,cfg_simulations, 0);
        
        MIs_chance_levels(isubject, ipoc) = nanmean(MIs_surr_distr{isubject, ipoc});
        POS_chance_levels(isubject, ipoc) = mean(POS_surr_distr{isubject, ipoc});
        z_rayleigh_chance_levels(isubject, ipoc) = mean(z_rayleigh_surr_distr{isubject, ipoc});
        U2watson_chance_levels(isubject, ipoc) = mean(U2watson_surr_distr{isubject, ipoc});
        rms_logregress_chance_levels(isubject, ipoc) = mean(rms_logregress_surr_distr{isubject, ipoc});
        
    end
end

% For the across-subject statistics later, compute distributions of
% surrogate average values.
MIs_surr_avg_distr = cell(1, length(cfg_simulations.strengths_phase_outcome_coupling));
POS_surr_avg_distr = cell(1, length(cfg_simulations.strengths_phase_outcome_coupling));
z_rayleigh_avg_surr_distr = cell(1, length(cfg_simulations.strengths_phase_outcome_coupling));
U2watson_surr_avg_distr = cell(1, length(cfg_simulations.strengths_phase_outcome_coupling));
rms_logregress_surr_avg_distr = cell(1, length(cfg_simulations.strengths_phase_outcome_coupling));
for inoise=1:length(cfg_simulations.strengths_phase_outcome_coupling)
    
    MIs_surr_avg_distr{inoise} = nan(1, cfg_simulations.nsurravg);
    POS_surr_avg_distr{inoise} = nan(1, cfg_simulations.nsurravg);
    z_rayleigh_avg_surr_distr{inoise} = nan(1, cfg_simulations.nsurravg);
    U2watson_surr_avg_distr{inoise} = nan(1, cfg_simulations.nsurravg);
    rms_logregress_surr_avg_distr{inoise} = nan(1, cfg_simulations.nsurravg);
    
    for isurr=1:cfg_simulations.nsurravg
        
        random_MI = nan(1, nsubjects);
        random_POS = nan(1, nsubjects);
        random_rayleigh = nan(1, nsubjects);
        random_U2 = nan(1, nsubjects);
        random_logregress = nan(1, nsubjects);
        for isubject=1:nsubjects

            random_MI(isubject) = datasample(MIs_surr_distr{isubject, inoise}, 1);
            random_POS(isubject) = datasample(POS_surr_distr{isubject, inoise}, 1);
            random_rayleigh(isubject) = datasample(z_rayleigh_surr_distr{isubject, inoise}, 1);
            random_U2(isubject) = datasample(U2watson_surr_distr{isubject, inoise}, 1);
            random_logregress(isubject) = datasample(rms_logregress_surr_distr{isubject, inoise}, 1);
            
        end
        MIs_surr_avg_distr{inoise}(isurr) = mean(random_MI);
        POS_surr_avg_distr{inoise}(isurr) = mean(random_POS);
        z_rayleigh_avg_surr_distr{inoise}(isurr) = mean(random_rayleigh);
        U2watson_surr_avg_distr{inoise}(isurr) = mean(random_U2);
        rms_logregress_surr_avg_distr{inoise}(isurr) = mean(random_logregress);
    end
end

if ~isdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_experiment' num2str(iexp) '.mat']), 'MIs_empirical');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_experiment' num2str(iexp) '.mat']), 'POS_empirical');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_empirical_experiment' num2str(iexp) '.mat']), 'z_rayleigh_empirical');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_experiment' num2str(iexp) '.mat']), 'U2watson_empirical');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical');

if ~isdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_experiment' num2str(iexp) '.mat']), 'POS_chance_levels');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_chance_levels_experiment' num2str(iexp) '.mat']), 'z_rayleigh_chance_levels');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels');

if ~isdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_experiment' num2str(iexp) '.mat']), 'POS_surr_distr');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_surr_distr_experiment' num2str(iexp) '.mat']), 'z_rayleigh_surr_distr');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr');

if ~isdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues_direct' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues_direct' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues_direct' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rayleigh_pvalues_direct_experiment' num2str(iexp) '.mat']), 'rayleigh_pvalues_direct');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues_direct' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_watson_pvalues_direct_experiment' num2str(iexp) '.mat']), 'watson_pvalues_direct');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues_direct' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_logregress_pvalues_direct_experiment' num2str(iexp) '.mat']), 'logregress_pvalues_direct');

if ~isdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep]))
end
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_avg_distr_experiment' num2str(iexp) '.mat']), 'MIs_surr_avg_distr');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_avg_distr_experiment' num2str(iexp) '.mat']), 'POS_surr_avg_distr');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_surr_avg_distr_experiment' num2str(iexp) '.mat']), 'z_rayleigh_avg_surr_distr');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_avg_distr_experiment' num2str(iexp) '.mat']), 'U2watson_surr_avg_distr');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_avg_distr_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_avg_distr');

if iexp==1     % save config-file only once
    save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep cfg_simulations.oscillation '_cfg_simulations.mat']), 'cfg_simulations');
end
end