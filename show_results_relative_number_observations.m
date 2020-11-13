function show_results_relative_number_observations(cfg_simulations, root_dir)
%SHOW_RESULTS_RELATIVE_NUMBER_OBSERVATIONS
%   Visualize sensitivity and False positive rate for the five statistical
%   tests as a function of relative number of observations for the two
%   groups.
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

% Load results of computations
% With resampling, sensitivity
load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_MI_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_MI_withresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_POS_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_POS_withresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_watson_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_watson_withresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_logregress_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_logregress_withresampling');

% Without resampling, sensitivity
load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_MI_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_MI_withoutresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_POS_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_POS_withoutresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_watson_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_watson_withoutresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_sensitivity_logregress_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'sensitivity_logregress_withoutresampling');

% With resampling, False Positive rate
load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_MI_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_MI_withresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_POS_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_POS_withresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_watson_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_watson_withresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'With_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_logregress_withresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_logregress_withresampling');

% Without resampling, False Positive rate
load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_MI_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_MI_withoutresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_POS_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_POS_withoutresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_watson_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_watson_withoutresampling');
load(strcat([root_dir 'Relative_trial_number' filesep 'Without_resampling' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_coupling_mode' num2str(cfg_simulations.coupling_mode) '_FPrate_logregress_withoutresampling_phaseoutcome_coupling' num2str(cfg_simulations.poc_relative_trial_number) '.mat']), 'FPrate_logregress_withoutresampling');

close all;

colors = {'r',[1 0.65 0],[0 0.5 0],'k',[0.12 0.56 1]};

figure('units','normalized','outerposition',[0.1645 0.2132 0.5605 0.5840]); set(gcf,'color','w');
l1=plot(FPrate_logregress_withresampling, '--', 'Color', colors{1}, 'Linewidth', 3); hold on;
p1=plot(FPrate_POS_withresampling, '--', 'Color', colors{2}, 'Linewidth', 3); hold on;
w1=plot(FPrate_watson_withresampling, '--', 'Color', colors{3}, 'Linewidth', 3); hold on;
m1=plot(FPrate_MI_withresampling, '--', 'Color', colors{5}, 'Linewidth', 3); hold on;
hold on;
l2=plot(FPrate_logregress_withoutresampling, 'Color', colors{1}, 'Linewidth', 3); hold on;
p2=plot(FPrate_POS_withoutresampling, 'Color', colors{2}, 'Linewidth', 3); hold on;
w2=plot(FPrate_watson_withoutresampling, 'Color', colors{3}, 'Linewidth', 3); hold on;
m2=plot(FPrate_MI_withoutresampling, 'Color', colors{5}, 'Linewidth', 3); hold on;
lgd=legend([l1 p1 w1 m1 l2 p2 w2 r2 m2], 'circ. log. regr. with resampling',  'POS with resampling', 'Watson with resampling', 'MI with resampling', 'circ. log. regr. without resampling', 'POS without resampling', 'Watson without resampling', 'MI without resampling', 'Location','northeast');
lgd.FontSize = 15;
hline(0.05, '--k');
ylim([0 0.1]);
xticks(1:length(cfg_simulations.proportion_group_sizes));
xticklabels({'20:80','30:70','40:60','50:50','60:40','70:30','80:20'});
yticks(0.01:0.01:0.1)
yticklabels({'1','2','3','4','5','6','7','8','9','10'});
ax = gca;
ax.FontSize = 16;
xlabel('Relative number of observations', 'FontSize', 20);
ylabel('% False Positives', 'FontSize', 20);
title('FP-rate', 'FontSize', 20);
grid on;

% sensitivity plot
figure('units','normalized','outerposition',[0.1645 0.2132 0.5605 0.5840]); set(gcf,'color','w');
l1=plot(sensitivity_logregress_withresampling, '--', 'Color', colors{1}, 'Linewidth', 3); hold on;
p1=plot(sensitivity_POS_withresampling, '--', 'Color', colors{2}, 'Linewidth', 3); hold on;
w1=plot(sensitivity_watson_withresampling, '--', 'Color', colors{3}, 'Linewidth', 3); hold on;
m1=plot(sensitivity_MI_withresampling, '--', 'Color', colors{5}, 'Linewidth', 3); hold on;
xticks(1:length(cfg_simulations.proportion_group_sizes));
xticklabels({'20:80','30:70','40:60','50:50','60:40','70:30','80:20'});
yticks(0.1:0.1:1)
yticklabels({'10','20','30','40','50','60','70','80','90','100'});
hold on;
l2=plot(sensitivity_logregress_withoutresampling, 'Color', colors{1}, 'Linewidth', 3); hold on;
p2=plot(sensitivity_POS_withoutresampling, 'Color', colors{2}, 'Linewidth', 3); hold on;
w2=plot(sensitivity_watson_withoutresampling, 'Color', colors{3}, 'Linewidth', 3); hold on;
m2=plot(sensitivity_MI_withoutresampling, 'Color', colors{5}, 'Linewidth', 3); hold on;
grid on;
lgd=legend([l1 p1 w1 m1 l2 p2 w2 m2], 'circ. log. regr. with resampling',  'POS with resampling', 'Watson with resampling', 'MI with resampling', 'circ. log. regr. without resampling', 'POS without resampling', 'Watson without resampling', 'MI without resampling', 'Location','northeast');
lgd.FontSize = 15;
ax = gca;
ax.FontSize = 16;
xlabel('Relative number of observations', 'FontSize', 20);
ylabel('% True Positives', 'FontSize', 20);
ylim([0 1]);
title(['Sensitivity, ' oscillation ', coupling mode ' num2str(cfg_simulations.coupling_mode) ', phase-outcome coupling = ' num2str(cfg_simulations.poc_relative_trial_number) ', ' num2str(cfg_simulations.nTrials_total_relative_trial_number) ' trials'], 'FontSize', 20);

end
