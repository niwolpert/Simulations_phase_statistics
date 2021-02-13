function [ POS, ITC_all, ITC_A, ITC_B ] = calc_POS( phases_hits, phases_misses )
%CALC_POS
%   Computes Phase Opposition Sum (POS).
%   Adapted from VanRullenR(2016) How to Evaluate Phase Differences between 
%   Trial Groups in Ongoing Electrophysiological Signals. Front. 
%   Neurosci.10:426. doi: 10.3389/fnins.2016.00426
%   
%   INPUTS
%   - phases_hits:          Phases of outcome condition 1 (in radians)
%   - phases_misses:        Phases of outcome condition 2 (in radians)
%
%   OUTPUTS
%   - POS:             Surrogate MI distribution
%   - ITC_all:             Surrogate POS distribution
%   - ITC_A:       Surrogate Rayleigh's z distribution
%   - ITC_B:        Surrogate Watson's U2 distribution

N = length(size(phases_hits));
if length(size(phases_misses))~=N
    fprintf('Exiting - phases_hits and phases_misses should have the same number of dimensions\n');
    return;
end;
allphases = cat(N,phases_hits,phases_misses);

%reformat data if necessary
if isreal(allphases(:)) %not complex numbers
%     fprintf('Not complex numbers, reformatting angles to complex values\n');
    allphases = exp(sqrt(-1).*allphases);
    phases_hits = exp(sqrt(-1).*phases_hits);
    phases_misses = exp(sqrt(-1).*phases_misses);
end;
if max(abs(abs(allphases(:))-1))>0.00001 %non-unit norm
    fprintf('Complex numbers with non-unit norm, normalizing\n');
    allphases = allphases ./ abs(allphases);
    phases_hits = phases_hits ./ abs(phases_hits);
    phases_misses = phases_misses ./ abs(phases_misses);
end;

ITC_all = squeeze(abs(mean(allphases,N)));
ITC_A = abs(mean(phases_hits,N));
ITC_B = abs(mean(phases_misses,N));
POS = squeeze(ITC_A+ITC_B)-2*ITC_all;

end

