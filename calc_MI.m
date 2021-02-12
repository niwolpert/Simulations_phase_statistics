function MI = calc_MI(phases_hits, phases_misses, nphasebins)
%CALC_PHMI
%   Calculates Modulation Index by Tort (2010).
%   
%   INPUTS
%   - phases_hits:          Phases of outcome condition 1 (in radians)
%   - phases_misses:        Phases of outcome condition 2 (in radians)
%   - nphasebins:           The number of bins into which phase is
%                           partitioned
%
%   OUTPUTS
%   - MI:                   Raw MI value
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

% extract number of hits and misses for each phase bin
[nhits_bins, ~] = histcounts(radtodeg(phases_hits), nphasebins);
[nmisses_bins, ~] = histcounts(radtodeg(phases_misses), nphasebins);
ntotal_bins=nhits_bins+nmisses_bins;

% compute hit rate for each phase bin
hr_bins = nhits_bins./ntotal_bins;

% If MI does not contain a hit in one bin (i.e. hit rate = 0), ignore this
% bin
if any(hr_bins==0)
    warning('Hit rate in one bin is zero, ignoring bin');
    
    hr_bins(hr_bins==0) = [];
    nphasebins = length(hr_bins);
end

MI = (log(nphasebins)-(-sum((hr_bins / sum(hr_bins)).* log((hr_bins /  sum(hr_bins))))))/log(nphasebins);

end

