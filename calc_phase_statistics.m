function [MI, POS, z_rayleigh, U2_watson, rms_logregress, p_rayleigh_direct, p_watson_direct, p_logregress_direct] = calc_phase_statistics(phases_hits, phases_misses, cfg_simulations, nresamples)
%CALC_PHASE_STATISTICS
%   Computes phase-outcome statistics using five circular tests: Modulation
%   Index (MI), Phase Opposition Sum (POS), Rayleigh test, Watson test and
%   Circular Logistic Regression. For the Rayleigh test, Watson test and
%   Circular Logistic Regression, additionally an individual p-value is
%   computed directly.
%   Optionally, applies a resampling procedure to control for imbalances in
%   the relative number of observations in the two conditions. Note that no
%   resampling is applied for the Rayleigh test as this is a one-sample
%   test.
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
%   - MI:                   Raw MI value
%   - POS:                  Raw POS value
%   - z_rayleigh:           z-value from Rayleigh-test
%   - U2_watson:            Watson's U2
%   - rms_logregress:       Root-mean-square of circular logistic
%                           regression
%   - p_rayleigh_direct:    Direct p-value from the Rayleigh test
%   - p_watson_direct:      Direct p-value from the Watson test
%   - p_logregress_direct:  Direct p-value from the circular logistic
%                           regression
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

% If resampling is used
if nresamples>0
    
    % initialize arrays for resample values
    all_POS_resample = nan(1, nresamples);
    all_MIs_resample = nan(1, nresamples);
    all_U2watson_resample = nan(1, nresamples);
    all_rms_logregress_resample = nan(1, nresamples);
    
    % identify condition with more trials
    ntrials = [length(phases_hits) length(phases_misses)];
    [~, group_max] = max(ntrials);
    ntrials_min = min(ntrials);
    
    % Case 1: more trials for hits than misses
    if group_max == 1

        for iresample=1:nresamples
            
            % sample as many trials from condition with more trials as 
            % there are trials in the condition with fewer trials
            phases_hits_resample = datasample(phases_hits, ntrials_min,'Replace',false);

            % compute POS on resampled set of hits
            all_POS_resample(iresample) = calc_POS(phases_hits_resample, phases_misses);
            
            % compute MI on resampled set of hits
            all_MIs_resample(iresample) = calc_MI(phases_hits_resample, phases_misses, cfg_simulations.nphasebins);
            
            % compute Watson's U2 on resampled set of hits
            all_U2watson_resample(iresample) = watsons_U2(phases_hits_resample', phases_misses');
            
            % compute root- mean square for circular regression on resampled set of hits
            all_rms_logregress_resample(iresample) = circlogress_calc_rms(phases_hits_resample,phases_misses);
            
        end

    % Case 2: more trials for misses than hits
    elseif group_max == 2
        
        for iresample=1:nresamples
            
            % sample as many trials from condition with more trials as 
            % there are trials in the condition with fewer trials
            phases_misses_resample  = datasample(phases_misses, ntrials_min,'Replace',false);
            
            % compute POS on resampled set of misses
            all_POS_resample(iresample) = calc_POS( phases_hits, phases_misses_resample );
            
            % compute MI on resampled set of misses
            all_MIs_resample(iresample) = calc_MI(phases_hits, phases_misses_resample, cfg_simulations.nphasebins);
            
            % compute Watson's U2 on resampled set of misses
            all_U2watson_resample(iresample) = watsons_U2(phases_hits', phases_misses_resample');
            
            % compute root- mean square for circular regression on resampled set of misses
            all_rms_logregress_resample(iresample) = circlogress_calc_rms(phases_hits,phases_misses_resample);
            
        end
    end
    
    % Estimate the test statistics as the mean of resample values
    POS = mean(all_POS_resample);
    MI = nanmean(all_MIs_resample);             % for MI, nans might occur due to missing trials in a phase bin
    U2_watson = mean(all_U2watson_resample);
    rms_logregress = mean(all_rms_logregress_resample); 
    
    % No resampling is applied for the Rayleigh test since it is computed
    % on one condition only
    [~, z_rayleigh] = circ_rtest(phases_hits);
    
% If no resampling is used
else
    MI = calc_MI(phases_hits, phases_misses, cfg_simulations.MI_nphasebins);
    [~, z_rayleigh] = circ_rtest(phases_hits);
    POS = calc_POS(phases_hits, phases_misses);
    U2_watson = watsons_U2(phases_hits', phases_misses');
    rms_logregress= circlogress_calc_rms(phases_hits,phases_misses);
end

% Additionally compute p-values directly for the Rayleigh test, Watson test 
% and logistic regression directly
p_rayleigh_direct = circ_rtest(phases_hits);
p_watson_direct = watsons_U2_approx_p(phases_hits', phases_misses');

% For circular logsitic regression, a p-value is computed by comparing the 
% full regression model with an intercept-only model
[~, pred_logregress] = circlogress_calc_rms(phases_hits,phases_misses);
outcomes = [zeros(length(phases_hits), 1); ones(length(phases_misses), 1)];
g1 = fitglm(pred_logregress,outcomes,'distribution','binomial','link','logit');
p_logregress_direct = coefTest(g1);

end
