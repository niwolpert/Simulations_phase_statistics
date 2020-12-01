%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SCRIPT_phase_stats_simulations_main
%%% 
%%% This script calls the functions needed to perform the simulations on
%%% circular statistics as reported in Wolpert & Tallon-Baudry, 2020.
%%% When using this function in any published study, please cite: Wolpert, 
%%% N. & Tallon-Baudry, C. (2020). Evaluation of different statistical 
%%% procedures to estimate coupling between oscillatory phase and 
%%% behavioral response (in press)
%%%
%%% The scripts and functions were written in Matlab version R2017b.
%%%
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%% USED TOOLBOXES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% In order for this script to run, make sure to install the following
%%% toolboxes:
%%% 1. Matlab's Parallel Computing Toolbox:
%%% https://fr.mathworks.com/products/parallel-computing.html
%%% 2. Circular statistics Toolbox: 
%%% https://fr.mathworks.com/matlabcentral/fileexchange/10676-circular-
%%% statistics-toolbox-directional-statistics
%%% Reference:
%%% Berens, P. (2009). CircStat: A MATLAB Toolbox for Circular Statistics. 
%%% J. Stat. Softw. 31, 1–21. doi:10.18637/jss.v031.i10
%%%
%%% Copyright (C) 2020, Laboratoire de Neurosciences Cognitives, Nicolai 
%%% Wolpert, Catherine Tallon-Baudry
%%% Email: nicolaiwolpert@gmail.com
%%% 
%%% DISCLAIMER:
%%% This code is provided without explicit or implicit guarantee, and 
%%% without any form of technical support. The code is not intended to be 
%%% used for clinical purposes. The functions are free to use and can be
%%% redistributed, modified and adapted, under the terms of the CC BY-NC-SA
%%% version of creative commons license (see
%%% <https://creativecommons.org/licenses/>).
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;

% Specify root directory containing scripts and functions
root_dir = '';

% Specify path to circular statistics toolbox
path_circstats = '';

addpath(genpath(root_dir));
addpath(genpath(path_circstats));

%% Specify parameters for simulations & statistics

% note them in a config structure
cfg_simulations = [];

% note sample rate of the oscillatory time series
cfg_simulations.fsample = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION-PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose if alpha (8-12 Hz) or EGG (~0.05 Hz) oscillation should be used
cfg_simulations.oscillation = 'alpha';

% Choose the "coupling mode" which specifies the simulated underlying 
% relationship between oscillatory phase and outcome, for which statistics
% will be estimated.
% 1: 1:1 coupling (events of the two outcomes each cluster at one given
% portion of the oscillatory cycle
% 2: 2:1 coupling (events of the two outcomes each cluster at two portions 
% of the oscillatory cycle
% 3: 3:1 coupling (events of the two outcomes each cluster at three 
% portions of the oscillatory cycle
% 4: 4:1 coupling (events of the two outcomes each cluster at four portions 
% of the oscillatory cycle
cfg_simulations.coupling_mode = 1;

% specify the range of the phase for defining probability functions
cfg_simulations.phase_range = -pi:0.00001:pi; 

% Create probabilities that specify the probabiliy of an outcome as a 
% function of oscillatory phase
if cfg_simulations.coupling_mode==1

    %% Mode 1: 1:1 coupling pattern 
    % Hits+misses clustered around one different phase each
    
    x = linspace(pi, pi*3, length(cfg_simulations.phase_range));
    cfg_simulations.p_hit = rescale(cos(x), 0.2, 0.8);
    cfg_simulations.p_miss = 1-cfg_simulations.p_hit;
    
elseif cfg_simulations.coupling_mode==2

    %% Mode 2: 2:1 coupling pattern
    % Hits and Misses clustered around two different phases each
    
    x = linspace(0, pi*4, length(cfg_simulations.phase_range));
    cfg_simulations.p_hit = cos(x);
    cfg_simulations.p_hit = rescale(cfg_simulations.p_hit, 0.2, 0.8);
    
    cfg_simulations.p_miss = 1-cfg_simulations.p_hit;

elseif cfg_simulations.coupling_mode==3
    
    %% Mode 3: 3:1 coupling pattern
    
    x = linspace(0, pi*6, length(cfg_simulations.phase_range));
    cfg_simulations.p_hit = cos(x);
    cfg_simulations.p_hit = rescale(cfg_simulations.p_hit, 0.2, 0.8);
    
    cfg_simulations.p_miss = 1-cfg_simulations.p_hit;
    
elseif cfg_simulations.coupling_mode==4
    
    %% Mode 4: 4:1 coupling pattern
    
    x = linspace(0, pi*8, length(cfg_simulations.phase_range));
    cfg_simulations.p_hit = cos(x);
    cfg_simulations.p_hit = rescale(cfg_simulations.p_hit, 0.2, 0.8);
    
    cfg_simulations.p_miss = 1-cfg_simulations.p_hit;
    
end
close all; fig_distr = figure('Color', 'w');
hits = plot(cfg_simulations.phase_range, cfg_simulations.p_hit, 'b', 'LineWidth', 2);
hold on;
misses = plot(cfg_simulations.phase_range, cfg_simulations.p_miss, 'r', 'LineWidth', 2);
title('Probability functions for outcomes', 'FontSize', 15);
xlabel('Phase (rad)', 'FontSize', 12);
ylabel('Outcome probability', 'FontSize', 12);
legend([hits misses], {'hits', 'misses'});
ax = gca;
ax.FontSize = 16;
xlim([-pi pi]); ylim([0 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% STATISTICAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%

% Select number of permutations to create surrogate distributions for the
% estimation of chance level and individual p-values. 100 permutations are
% generally sufficient.
cfg_simulations.nperm = 100;        

% Number of reshuffles to create a distribution of surrogate averages, for
% estimating significance on the group level.
cfg_simulations.nsurravg = 1000;

% MI parameters:
% number of bins into which EGG phase is partitioned
cfg_simulations.MI_nphasebins = 10;
% edges of the bins
cfg_simulations.bin_edges = -pi:(2*pi)/cfg_simulations.MI_nphasebins:pi;   

%% Load oscillatory phase time series

if strcmp(cfg_simulations.oscillation, 'alpha')
    
    fprintf('\n##### Loading alpha phase time series...\n');
    
    phase_all_subjects = cell(1, 30);
    
    for isubject=1:30
        
        load(strcat([root_dir cfg_simulations.oscillation filesep 'Subject' num2str(isubject,'%02.f') '_alpha_phase.mat']));
        phase_all_subjects{isubject} = alpha_phase;
        
    end
    
elseif strcmp(cfg_simulations.oscillation, 'EGG')
    
    fprintf('\n##### Loading EGG phase time series...\n');
    
    phase_all_subjects = cell(1, 30);
    
    for isubject=1:30
        
        load(strcat([root_dir cfg_simulations.oscillation filesep 'Subject' num2str(isubject,'%02.f') '_EGG_phase.mat']));
        phase_all_subjects{isubject} = alpha_phase;
        
    end
end

%% Compare statistical tests & methods to assess significance at the group level
% Compute statistics for the five statistical tests over different
% strengths of phase outcome coupling. The results of these computations
% will later be used for the following analyses:
% 1. Sensitivity of different tests across different coupling modes
% 2. Different ways of assessing significance at the group level
% 3. Estimating chance level as the mean vs. median of surrogate
%    distributions
% 4. Comparing permutation statistics with tabulated statistics

% Define the strengths of phase-outcome coupling for which statistics will 
% be estimated
% E.g. for a phase-outcome coupling of 0.2, outcomes for 20% of trials will 
% be assigned as a function of phase while for the remaining trials,
% outcome will be randomly assigned
cfg_simulations.strengths_phase_outcome_coupling = 0:5:50;

% Fix number of trials
cfg_simulations.nTrials_poc = 250;

% specify ISI range
ISI_range=[0.1 2];

% choose number of virtual experiments
nexperiments=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Launch computations %%%%%%%%%%%%%%%%%%%%%%%%%%%
% In a first step, phase-outcome statistics (POS, MI, ...) are computed for
% each phase-outcome coupling strength

% open parallel pool
parpool local

parfor iexp=1:nexperiments
    
    comp_stats_per_poc(cfg_simulations, iexp, phase_all_subjects, root_dir, ISI_range);
    
end

% close parallel pool
poolobj = gcp('nocreate');
delete(poolobj);

%%% Compute sensitivity and False Positive rate based on results
% Results are saved in root directory
signific_thresh = 0.05;
merge_results_per_phaseoutcome_coupling(cfg_simulations, root_dir, nexperiments, signific_thresh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Evaluate results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% 1. Sensitivity of different tests across different coupling modes

% Choose method to estimate significance on the group level
% Options: 'empirical_vs_chance', 'surravg', 'combined_stouffer', 
% 'combined_fisher' or 'combined_edington'
method_estimate_significance = 'empirical_vs_chance';
show_results_comparison_tests(cfg_simulations, root_dir, method_estimate_significance)

%%% 2. Compare methods to estimate signfificance on the group-level

show_results_comparison_methods_grouplevelstats(cfg_simulations, root_dir)

%%% 3. Estimate chance level as the mean vs. median of surrogate 
%%% distributions

show_results_chance_level_mean_vs_median(cfg_simulations, root_dir, nexperiments)

%%% 4. Compare permutation statistics with tabulated statistics

show_results_comparison_permutation_vs_tabulated_stats(cfg_simulations, root_dir)

%% Number of trials in total
% Systematically increment the number of trials and estimate sensitivity
% and False Positive rate for each statistical test

% Specify levels of trial numbers
cfg_simulations.levels_nTrials = [50 100 150 200 250 300 350 400 450 500];

% Fix strength of phase-outcome coupling
cfg_simulations.poc_simulations_Ntrials = 20;

% Specify ISI range
% Since many trials are included in the last steps, we need to make sure
% that enough 'events' are created by narrowing a bit the ISI distribution
ISI_range=[0.1 1.5];

% choose number of virtual experiments
nexperiments=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Launch computations %%%%%%%%%%%%%%%%%%%%%%%%%%%
% In a first step, phase-outcome statistics (POS, MI, ...) are computed for
% each number of trials

% open parallel pool
parpool local

parfor iexp=1:nexperiments
    
    comp_stats_per_nTrials(cfg_simulations, iexp, phase_all_subjects, root_dir, ISI_range);
    
end

% close parallel pool
poolobj = gcp('nocreate');
delete(poolobj);

%%% Compute sensitivity and False Positive rate based on results
% Results are saved in root directory
signific_thresh = 0.05;
merge_results_per_ntrials(cfg_simulations, root_dir, nexperiments, signific_thresh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Evaluate results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

show_results_number_of_trials(cfg_simulations, root_dir)

%% Relative number of observations
% Systematically vary the relative number of observations for hits and
% misses while keeping the total number of trials and phase-outcome 
% coupling strength constant. Use a resampling vs. no resampling procedure
% to control for trial imbalance to investigate how it affects sensitivity
% and False Positive rate.

% Choose levels of relative number of observations (e.g. '0.2' = 20% hits,
% 80% misses)
cfg_simulations.proportion_group_sizes = [0.2 0.3 0.4 0.5 0.6 0.7 0.8];

% Specify number of iterations for resampling procedure. 
cfg_simulations.nresamples_relative_trial_number = 100;

% Fix phase-outcome coupling strength
cfg_simulations.poc_relative_trial_number = 15;

% Fix total number of trials
cfg_simulations.nTrials_total_relative_trial_number = 300;

% Specify ISI range
% Since many trials are included in the last steps, we need to make sure
% that enough 'events' are created by narrowing a bit the ISI distribution
ISI_range = [0.1 1.5];

% choose number of virtual experiments
nexperiments=1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Launch computations %%%%%%%%%%%%%%%%%%%%%%%%%%%
% In a first step, phase-outcome statistics (POS, MI, ...) are computed for
% each relative number of observations

% open parallel pool
parpool local

parfor iexp=200:nexperiments
    
    comp_stats_per_relative_number_observations(cfg_simulations, iexp, phase_all_subjects, root_dir, ISI_range);

end

% close parallel pool
poolobj = gcp('nocreate');
delete(poolobj);

%%% Compute sensitivity and False Positive rate based on results
% Results are saved in root directory
thresh_sign = 0.05;
merge_results_per_relative_number_observations(cfg_simulations, root_dir, nexperiments, thresh_sign)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Evaluate results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

show_results_relative_number_observations(cfg_simulations, root_dir);

%% Optimal number of bins for Modulation Index (MI)
% Because MI depends on the logarithm of hit rate, it cannot be computed in
% the case that a given phase bin contains no hit. Here, estimate the 
% relationship between the total number of trials (hits+misses), number of 
% phase bins and the likelihood that MI cannot be computed.
% Systematically vary the number of bins for MI and the total number of
% trials in 1,000 virtual "subjects", and compute the propotion of subjects
% in which MI could not be computed.

% Fix phase-outcome coupling strength (not important for this simulation)
phase_outcome_coupling = 30;

% Specify levels of bin numbers for MI
levels_nbins = [5 8 10 13 17 20];

% Specify levels of trials in total
levels_nTrials = [50 75 100 125 150 175 200 225 250 275 300 325 350 400 450 500];

ISI_range=[0.1 1.5];

% Select number of virtual subjects
nsubjects=1000;

% This matrix specifies if there is any bin with no hit (1 if this is
% the case, else 0), for each subject (1st dimension), number of trials
% (2nd dimension and number of bins (3rd dimension)
empty_bin_per_subject_nbins_trialnumber = zeros(nsubjects, length(levels_nTrials), length(levels_nbins));

for isubject=1:nsubjects
    
    fprintf(['\n### Subject' num2str(isubject) '...\n']);
    
    % Here, with each virtual subject, we use a random phase time series of
    % the 30 phase time series from real participants.
    % Start out with the maximum number of trials and later create
    % subsamples to systematically decrease the number of trials
    [phases_hits_original, phases_misses_original] = distribute_outcomes(cfg_simulations, ISI_range, phase_all_subjects{randi(length(phase_all_subjects))}, phase_outcome_coupling, max(levels_nTrials));
    % mix up phases
    phases_hits_original = datasample(phases_hits_original, length(phases_hits_original), 'Replace', false);
    phases_misses_original = datasample(phases_misses_original, length(phases_misses_original), 'Replace', false);
    
    for iTrials=1:length(levels_nTrials)
        
        nTrials = levels_nTrials(iTrials);
        
        % Subsample a set of trials for hits and misses to have exactly as
        % many trials as specificed
        phases_hits = phases_hits_original(1:round(nTrials/2));
        phases_misses = phases_misses_original(1:round(nTrials/2));
        
        for ibin=1:length(levels_nbins)
            
            nphasebins = levels_nbins(ibin);
            
            % If MI cannot be computed due to a missing hit in at least one
            % phase bin, it results in an nan
            if isnan(calc_MI(phases_hits, phases_misses, nphasebins))
                
                empty_bin_per_subject_nbins_trialnumber(isubject, iTrials, ibin) = 1;
                
            end
        end
    end
end

% Save intermediate results
if ~isdir(strcat([root_dir 'Optimal_number_bins_MI' filesep]))
    mkdir(strcat(strcat([root_dir 'Optimal_number_bins_MI' filesep])))
end
save(strcat([root_dir 'Optimal_number_bins_MI' filesep 'empty_bin_per_subject_nbins_trialnumber_' cfg_simulations.oscillation '_oscillation.mat']), 'empty_bin_per_subject_nbins_trialnumber');

% Compute proportion of subjects that would have to be discarded per trial
% number and bin number
perc_missing_per_bin = flipud(rot90(squeeze(sum(empty_bin_per_subject_nbins_trialnumber, 1))/nsubjects))*100;

% Save final results
save(strcat([root_dir 'Optimal_number_bins_MI' filesep 'perc_missing_per_bin_' cfg_simulations.oscillation '_oscillation.mat']), 'perc_missing_per_bin');

%%% Show results: percentage of subjects where MI cannot be computed, as a
%%% function of the number of trials and phase bins
close all;
load(strcat([root_dir 'Optimal_number_bins_MI' filesep 'perc_missing_per_bin_' cfg_simulations.oscillation '_oscillation.mat']), 'perc_missing_per_bin');
figure('units','normalized','outerposition',[0.1105 0.2229 0.5590 0.6083]); set(gcf,'color','w');
colors = {[0 0 1], [0 0.5 1], [0 1 1], [0 1 0], [1 0.5 0], [1 0 0]};
x=(1:length(levels_nTrials));
l=nan(1, length(levels_nbins));
symbols = {'x','v','s','*','o','^'};
for ibin=1:length(levels_nbins)
    hold on;
    
    l(ibin) = plot(x, perc_missing_per_bin(ibin, :), 'LineWidth', 2, 'Color', colors{ibin});
    
    hold on;
    plot(x, perc_missing_per_bin(ibin, :), symbols{ibin}, 'LineWidth', 1.5, 'MarkerSize', 10, 'Color', colors{ibin});
    
end
ax = gca;
ax.FontSize = 16;
% yticks(0:0.1:1);
xticks(1:length(levels_nTrials));
xticklabels(levels_nTrials);
grid on;
xlim([0.5 length(levels_nTrials)+0.3]);
lgd=legend(fliplr(l), fliplr({'5 bins','8 bins','10 bins','13 bins','17 bins','20 bins'}));
lgd.FontSize = 15;
xlabel('#trials', 'FontSize', 20);
ylabel('%participants with nan', 'FontSize', 20); xlim([0.5 13.5]);
title('Percentage of subjects where MI cannot be computed', 'FontSize', 20);
