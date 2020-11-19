function show_results_number_of_trials(cfg_simulations, root_dir)
%SHOW_RESULTS_NUMBER_OF_TRIALS
%   Shows the results of simulations on the number of trials.
%   
%   INPUTS
%   - cfg_simulations:                 Configuration structure with
%                                      simulation parameters
%   - root_dir:                        Root directory where results will be
%                                      saved
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

load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_FPrate_ntrials_MI_phase_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'FPrate_ntrials_MI');
load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_FPrate_ntrials_POS_noise_level' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'FPrate_ntrials_POS');
load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_FPrate_ntrials_z_rayleigh_noise_level' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'FPrate_ntrials_rayleigh');
load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_FPrate_ntrials_U2watson_noise_level' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'FPrate_ntrials_watson');
load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_FPrate_ntrials_rms_logregress_noise_level' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'FPrate_ntrials_logregress');

load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_sensitivity_ntrials_MI_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'sensitivity_ntrials_MI');
load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_sensitivity_ntrials_POS_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'sensitivity_ntrials_POS');
load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_sensitivity_ntrials_z_rayleigh_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'sensitivity_ntrials_rayleigh');
load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_sensitivity_ntrials_U2watson_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'sensitivity_ntrials_watson');
load(strcat([root_dir 'Number_of_trials' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_sensitivity_ntrials_rms_logregress_outcome_coupling' num2str(cfg_simulations.poc_simulations_Ntrials) '.mat']), 'sensitivity_ntrials_logregress');

close all;

colors = {'r',[1 0.65 0],[0 0.5 0],'k',[0.12 0.56 1]};

figure('units','normalized','outerposition',[0.1488 0.1035 0.6148 0.7840]); set(gcf,'color','w');
subplot(2,1,1);
l=plot(sensitivity_ntrials_logregress, 'Color', colors{1}, 'Linewidth', 3); hold on;
p=plot(sensitivity_ntrials_POS, 'Color', colors{2}, 'Linewidth', 3); hold on;
w=plot(sensitivity_ntrials_watson, 'Color', colors{3}, 'Linewidth', 3); hold on;
r=plot(sensitivity_ntrials_rayleigh, 'Color', colors{4}, 'Linewidth', 3); hold on;
m=plot(sensitivity_ntrials_MI, 'Color', colors{5}, 'Linewidth', 3); hold on;
xticks(1:length(cfg_simulations.levels_nTrials));
xticklabels({cfg_simulations.levels_nTrials});
lgd=legend([l p w r m], 'Circ. log. regr.', 'POS', 'Watson', 'Rayleigh', 'MI', 'Location','northeast');
yticks(0:0.1:1);
lgd.FontSize = 15;
ax = gca;
ax.FontSize = 16;
xlim([1 8]);
ylim([0 1]);
xlabel('Number of trials', 'FontSize', 20);
ylabel('% True Positives', 'FontSize', 20);
title(['Sensitivity, ' cfg_simulations.oscillation ', coupling mode ' num2str(cfg_simulations.coupling_mode) ', phase-outcome coupling = ' num2str(cfg_simulations.poc_simulations_Ntrials) '%'], 'FontSize', 20);
% grid on;

subplot(2,1,2);
l=plot(FPrate_ntrials_logregress,  'Color', colors{1}, 'Linewidth', 3); hold on;
p=plot(FPrate_ntrials_POS, 'Color',colors{2}, 'Linewidth', 3); hold on;
w=plot(FPrate_ntrials_watson, 'Color', colors{3}, 'Linewidth', 3); hold on;
r=plot(FPrate_ntrials_rayleigh, 'Color', colors{4}, 'Linewidth', 3); hold on;
m=plot(FPrate_ntrials_MI, 'Color', colors{5}, 'Linewidth', 3); hold on;
xticks(1:length(cfg_simulations.levels_nTrials));
xticklabels({cfg_simulations.levels_nTrials});
ax = gca;
ax.FontSize = 16;
xlim([1 8]);
ylim([0 0.1]);
hline(0.05, '--k');
xlabel('Number of trials', 'FontSize', 20);
ylabel('% False Positives', 'FontSize', 20);
title(['FP-rate, ' cfg_simulations.oscillation ', coupling mode ' num2str(cfg_simulations.coupling_mode) ', phase-outcome coupling = ' num2str(cfg_simulations.poc_simulations_Ntrials) '%'], 'FontSize', 20);
% grid on;

end
