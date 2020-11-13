function [phases_hits,phases_misses,ind_hits,ind_misses,events] = distribute_outcomes(cfg_simulations, ISI_range, phase_time_series_individual, phase_outcome_coupling, nTrials_wanted)
%DISTRIBUTE_OUTCOMES
%   Creates a time series of 'stimuli', 'trials' or 'events' with two 
%   outcome conditions (here called 'hits' and 'misses'), based on 
%   parameters specified in 'cfg_simulations'. These parameters include the 
%   probability function for outcomes by phase, the number of trials, and 
%   the strength of the phase-outcome coupling.
%   In a first step, events are distributed with a random interval selected
%   from a flat distribution in the range specified by 'ISI_range'. Next,
%   hits and misses are assigned to these events based on the probability 
%   functions for outcomes by phase (in fields 'p_hit' and 'p_miss'). 
%   Finally, the phase-outcome coupling strength is controlled by selecting 
%   a percentage of trials (100-phase_outcome_coupling) where outcome is 
%   randomly reassigned, independently of phase.
%
%   INPUTS
%   - cfg_simulations:                 Configuration structure with
%                                      simulation parameters
%   - ISI_range:                       The lower and upper limit for
%                                      intervals between events
%   - phase_time_series_individual:    The phase time series of a given
%                                      subject
%   - phase_outcome_coupling:          The strength of the phase outcome
%                                      coupling (integer between 0 and 100)
%   - nTrials_wanted:                  Number of trials in total
%
%   OUTPUTS
%   - phases_hits:                     Phases for hits in radians
%   - phases_misses:                   Phases for misses in radians
%   - ind_hits:                        Sample indeces for hits with respect
%                                      to the individual phase time series
%   - ind_misses:                      Sample indeces for hits with respect
%                                      to the individual phase time series
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

fprintf('\nCreating time series of events of two outcome conditions...\n');

% Compute the percentage of trials where outcome becomes assigned
% independently of phase ("noise trials")
percentage_noise = 100-phase_outcome_coupling;

%%% Assign events
% this matrix stores events of outcomes 'hits' and 'misses, represented by 
% 1 and -1 respectively
events = zeros(1, length(phase_time_series_individual));

% Distribute some events by noting the samples of the respective event
% onsets, using a flat ISI distribution.
% In the beginning, this will results in a bit more trials than wanted
ISI_range_samples=ISI_range.*1000;
isample=randi(ISI_range_samples);
iTrial = 0;
ind_stim = [];      % stores all samples of artificial stimulus onsets
while isample<length(phase_time_series_individual)
    
    iTrial = iTrial+1;
    
    ind_stim = [ind_stim isample];
    
    isample = isample+randi(ISI_range_samples);
    
end
% Compute the actual number of trials obtained so far
nTrials_actual=length(ind_stim);
if nTrials_actual<nTrials_wanted
    error('To few trials obtained. Try changing the ISI distribution to obtain more trials');
end

% Rotate probability function for hits and misses such that each
% participant has a different preferred phase
shift = randi([1 length(cfg_simulations.phase_range)]);
cfg_simulations.p_hit = circshift(cfg_simulations.p_hit, shift);
cfg_simulations.p_miss = circshift(cfg_simulations.p_miss, shift);

% Assign hits and misses to each trial, with probabilities defined 
% according to the phase at which the event occurs.
for iTrial = 1:nTrials_actual
    
    % note phase of current trial
    phase = phase_time_series_individual(ind_stim(iTrial));
    
    % find index of phase in phase_range that closest matches the phase
    ind_phase = nearestpoint(phase, cfg_simulations.phase_range);
    
    % assign a hit or miss to this stimulus according to the 
    % probabilities for hits and misses for this given phase
    events(ind_stim(iTrial)) = prob_select([cfg_simulations.p_miss(ind_phase) cfg_simulations.p_hit(ind_phase)],[-1 1]);
    
end

% Compute the number of noise and no-noise trials
nNoiseTrials = round(nTrials_actual*(percentage_noise/100));
noise_trials = datasample(ind_stim, nNoiseTrials, 'Replace',false);

% Add noise:
% Make a proportion of trials, controlled by "noise_level", to be not
% dependent on phase but on noise. This means that for these trials, we
% assign events 1 and 2 randomly without taking phase into account.
for iNoiseTrial = 1:nNoiseTrials
    
    events(noise_trials(iNoiseTrial)) = prob_select([0.5 0.5],[-1 1]);

end

% Remove events to have exactly as many trials as specified in
% 'nTrials_wanted' and with equal amount of hits and misses
ind_hits = find(events==1);
ind_misses = find(events==-1);
nhits=length(ind_hits);
nmisses=length(ind_misses);

ind_events_to_remove = [datasample(ind_hits, nhits-(nTrials_wanted/2), 'Replace',false) datasample(ind_misses, nmisses-(nTrials_wanted/2), 'Replace',false)];
ind_stim = setdiff(ind_stim, ind_events_to_remove);

events(ind_events_to_remove)=0;

ind_hits = find(events==1);
ind_misses = find(events==-1);
nhits=length(ind_hits);
nmisses=length(ind_misses);

% note time series of hits and misses


% note phases of events 1 and 2
phases_hits = phase_time_series_individual(ind_hits);
phases_misses = phase_time_series_individual(ind_misses);

fprintf(['Created ' num2str(length(phases_hits)) ' hits and ' num2str(length(phases_misses)) ' misses with ' num2str(phase_outcome_coupling) '%% phase-outcome coupling strength.\n']);

end
