function [rms, predictors] = circlogress_calc_rms(phases_hits,phases_misses)
%CIRCLOGRESS_CALC_RMS
%   Calculates root mean square for circular logistic regression.
%   Adapted from Zoefel, B, Davis, M.H., Valente, G., Riecke, L. (2019).
%   How to test for phasic modulation of neural and behavioural responses.
%   NeuroImage, Volume 202, 116175, ISSN 1053-8119, 
%   https://doi.org/10.1016/j.neuroimage.2019.116175.
%
%   INPUTS
%   - phases_hits:          Phases of outcome condition 1 (in radians)
%   - phases_misses:        Phases of outcome condition 2 (in radians)
%
%   OUTPUTS
%   - rms:                  Root-mean-square of obtained predictor
%                           coefficients
%   - predictors:           Predictor coefficients

phases_all = [phases_hits phases_misses]';
outcomes = [ones(length(phases_hits), 1); repmat(2, length(phases_misses), 1)];
predictors = zeros(size(phases_all,1),2);
predictors(:,1) = sin(phases_all);
predictors(:,2) = cos(phases_all);
betas = mnrfit(predictors,outcomes);
rms = sqrt(betas(2)^2 + betas(3)^2);

end