function [ element_selected ] = prob_select(probs,elements)
%PROB_SELECT
%   Select one of two elements based on probabilities for each element.
%
%   INPUTS
%   - probs:                1xN array of probabilities for each of the N 
%                           elements (sums up to 1)
%   - elements:             1xN array of N elements among which to choose 
%                           an element
%
%   OUTPUTS
%   - element_selected:     The chosen element
%
% Nicolai Wolpert, 2020

rng('shuffle');

p = cumsum([0; probs(1:end-1).'; 1+1e3*eps]);
[a a] = histc(rand,p);
element_selected = elements(a); 
    
end