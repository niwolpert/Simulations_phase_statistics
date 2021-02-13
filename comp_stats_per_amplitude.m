function comp_stats_per_amplitude(cfg_simulations, iexp, pinknoise_all_subjects, sinewaves_all_subjects, root_dir, ISI_range)
%COMP_STATS_PER_AMPLITUDE
%   Runs simulations to compare the sensitivity for different circular 
%   statistical tests (Phase Opposition Sum (POS), Circular Logistic 
%   Regression, Watson test, Modulation Index (MI) and Rayleigh test) as a
%   function of the amplitude of the underlying oscillation that modulates
%   outcome.
%   First create a synthetic 10 Hz oscillation as a sinewave of amplitude 
%   scaled to [-1;1], with a sampling frequency of 1000 Hz and 15 minutes 
%   duration. Then, for each subject, assign hits and misses based on the 
%   synthetic 10Hz sinewave. Then modulate the amplitude of the 10Hz 
%   oscillation by a scaling factor ranging between 0 and 0.2 before adding
%   it to background noise, generated as pink noise with an amplitude 
%   rescaled to [-1, 1]. The resulting combined signal is then filtered 
%   around 10 Hz (±1) using a 6th order Butterworth zero-phase shift 
%   filter, and the Hilbert transform is applied on the combined signal to 
%   extract instantaneous phase. The resulting phase time series thus
%   represents the "empirical" phase time series whose signal-to-noise 
%   ratio depends on the amplitude of the true underlying oscillation. 
%   Phases for hits and misses are extracted, and the phase-outcome 
%   statistics computed for each amplitude. This is repeated in a number of 
%   virtual experiments, to compute sensitivity and False Positive rate.
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

nsubjects = length(pinknoise_all_subjects);

MIs_empirical = nan(nsubjects, length(cfg_simulations.scaling_factors));
POS_empirical = nan(nsubjects, length(cfg_simulations.scaling_factors));
U2watson_empirical = nan(nsubjects, length(cfg_simulations.scaling_factors));
rms_logregress_empirical = nan(nsubjects, length(cfg_simulations.scaling_factors));

MIs_surr_distr = cell(nsubjects, length(cfg_simulations.scaling_factors));
POS_surr_distr = cell(nsubjects, length(cfg_simulations.scaling_factors));
U2watson_surr_distr = cell(nsubjects, length(cfg_simulations.scaling_factors));
rms_logregress_surr_distr = cell(nsubjects, length(cfg_simulations.scaling_factors));

MIs_chance_levels = nan(nsubjects, length(cfg_simulations.scaling_factors));
POS_chance_levels = nan(nsubjects, length(cfg_simulations.scaling_factors));
U2watson_chance_levels = nan(nsubjects, length(cfg_simulations.scaling_factors));
rms_logregress_chance_levels = nan(nsubjects, length(cfg_simulations.scaling_factors));

for isubject=1:nsubjects

    fprintf(['\n### Subject' num2str(isubject) '\n']);
    
    % Take sine wave that represents true underlying oscillation
    sinewave = sinewaves_all_subjects{isubject};
    
    % Extract phase
    phase_real = angle(hilbert(sinewave'))';
    
    % Create trials with equal amount of hits and misses
    [~,~,ind_events1,ind_events2] = distribute_outcomes(cfg_simulations, ISI_range, phase_real, cfg_simulations.amplitude_poc, cfg_simulations.amplitude_nTrials);
    
    % Compute stats for each amplitude of the sine wave
    for iamp=1:length(cfg_simulations.scaling_factors)
        
        % Superimpose sine wave on pink noise with given amplitude
        noise_plus_sinewave = pinknoise_all_subjects{isubject}+sinewave*cfg_simulations.scaling_factors(iamp);
        
        % Filter
        noise_plus_sinewave_filtered = ft_preproc_bandpassfilter(noise_plus_sinewave,cfg_simulations.fsample,[8 12]);
        
        % Compute phase
        phase_noise_plus_sinewave = angle(hilbert(noise_plus_sinewave_filtered'))';
        
        phases_hits=phase_noise_plus_sinewave(ind_events1);
        phases_misses=phase_noise_plus_sinewave(ind_events2);
        
        [MIs_empirical(isubject, iamp), POS_empirical(isubject, iamp), U2watson_empirical(isubject, iamp), rms_logregress_empirical(isubject, iamp)] = calc_phase_statistics(phases_hits, phases_misses, cfg_simulations, 0);
        
        [MIs_surr_distr{isubject, iamp}, POS_surr_distr{isubject, iamp}, U2watson_surr_distr{isubject, iamp}, rms_logregress_surr_distr{isubject, iamp}] = calc_chance_level_statistics(phases_hits,phases_misses,cfg_simulations,0);
        
        MIs_chance_levels(isubject, iamp) = nanmean(MIs_surr_distr{isubject, iamp});
        POS_chance_levels(isubject, iamp) = mean(POS_surr_distr{isubject, iamp});
        U2watson_chance_levels(isubject, iamp) = mean(U2watson_surr_distr{isubject, iamp});
        rms_logregress_chance_levels(isubject, iamp) = mean(rms_logregress_surr_distr{isubject, iamp});
        
    end
end

if ~isdir(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep]))
   
    mkdir(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep]))
    
end

save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_MIs_empirical_iexp' num2str(iexp) '.mat']), 'MIs_empirical');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_POS_empirical_iexp' num2str(iexp) '.mat']), 'POS_empirical');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_U2watson_empirical_iexp' num2str(iexp) '.mat']), 'U2watson_empirical');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep 'amplitudes_rms_logregress_empirical_iexp' num2str(iexp) '.mat']), 'rms_logregress_empirical');

save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep 'amplitudes_MIs_chance_levels_iexp' num2str(iexp) '.mat']), 'MIs_chance_levels');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep 'amplitudes_POS_chance_levels_iexp' num2str(iexp) '.mat']), 'POS_chance_levels');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep 'amplitudes_U2watson_chance_levels_iexp' num2str(iexp) '.mat']), 'U2watson_chance_levels');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep 'amplitudes_rms_logregress_chance_levels_iexp' num2str(iexp) '.mat']), 'rms_logregress_chance_levels');

save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep 'amplitudes_MIs_surr_distr_iexp' num2str(iexp) '.mat']), 'MIs_surr_distr');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep 'amplitudes_POS_surr_distr_iexp' num2str(iexp) '.mat']), 'POS_surr_distr');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep 'amplitudes_U2watson_surr_distr_iexp' num2str(iexp) '.mat']), 'U2watson_surr_distr');
save(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep 'amplitudes_rms_logregress_surr_distr_iexp' num2str(iexp) '.mat']), 'rms_logregress_surr_distr');

end
