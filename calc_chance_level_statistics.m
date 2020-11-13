function [MIs_surr, POS_surr, zrayleigh_surr, U2watson_surr, rms_logregress_surr] = calc_chance_level_statistics(phases_hits, phases_misses, cfg_simulations, nresamples)
%CALC_CHANCE_LEVEL_STATISTICS
%   Computes chance level statistics for Modulation Index (MI), Phase
%   Opposition Sum (POS), Rayleigh test, Watson test and Circular Logistic
%   Regression, using a permutation-based approach. Distributions of phase
%   statistics under the null hypothesis of no phase-outcome relationship
%   are computed by randly reshuffling outcomes.
%   
%   INPUTS
%   - phases_hits:          Phases of outcome condition 1 (in radians)
%   - phases_misses:        Phases of outcome condition 2 (in radians)
%   - cfg_simulations:      Configuration structure with simulation 
%                           parameters
%   - nresamples:           Number of resamples for the procedure to
%                           control for imbalances in the relative number 
%                           of trials. (If zero, no resampling is applied)
%
%   OUTPUTS
%   - MIs_surr:             Surrogate MI distribution
%   - POS_surr:             Surrogate POS distribution
%   - zrayleigh_surr:       Surrogate Rayleigh's z distribution
%   - U2watson_surr:        Surrogate Watson's U2 distribution
%   - rms_logregress_surr:  Surrogate distribution of logistic regression
%                           RMS values
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

% Initialize arrays for surrogate distributions
POS_surr = nan(1, cfg_simulations.nperm);
MIs_surr = nan(1, cfg_simulations.nperm);
zrayleigh_surr = nan(1, cfg_simulations.nperm);
U2watson_surr = nan(1, cfg_simulations.nperm);
rms_logregress_surr = nan(1, cfg_simulations.nperm);

% Perform as many permutations as specified in 'cfg_simulations.nperm'
dispstat('','init');
for isurr=1:cfg_simulations.nperm
    
    dispstat(sprintf('Permutation %d%',isurr));
    
    % randomly permute labels
    N = length(size(phases_hits));

    allphases = cat(N,phases_hits,phases_misses);
    
    % create permutation indeces
    order = randperm(size(allphases,N));
    perm1 = order(1:size(phases_hits,N));
    perm2 = order(size(phases_hits,N)+1:end);
    
    phases_hits_perm = allphases(perm1);
    phases_misses_perm = allphases(perm2);
    
    [MIs_surr(isurr), POS_surr(isurr), zrayleigh_surr(isurr), U2watson_surr(isurr), rms_logregress_surr(isurr)] = calc_phase_statistics(phases_hits_perm, phases_misses_perm, cfg_simulations, nresamples);
    
end
end