function show_results_comparison_methods_grouplevelstats(cfg_simulations, root_dir)
%SHOW_RESULTS_COMPARISON_METHODS_GROUPLEVELSTATS
%   Shows the results on comparison of sensitivity and False positive rate 
%   of the five different methods to estimate significance on the group
%   level.
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

% FP rate

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_MI.mat']), 'fp_empirical_vs_chance_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_POS.mat']), 'fp_empirical_vs_chance_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_z_rayleigh.mat']), 'fp_empirical_vs_chance_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_U2watson.mat']), 'fp_empirical_vs_chance_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_rms_logregress.mat']), 'fp_empirical_vs_chance_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_MI.mat']), 'fp_surravg_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_POS.mat']), 'fp_surravg_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_z_rayleigh.mat']), 'fp_surravg_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_U2watson.mat']), 'fp_surravg_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_rms_logregress.mat']), 'fp_surravg_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_MI.mat']), 'fp_combined_stouffer_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_POS.mat']), 'fp_combined_stouffer_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_z_rayleigh.mat']), 'fp_combined_stouffer_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_U2watson.mat']), 'fp_combined_stouffer_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_rms_logregress.mat']), 'fp_combined_stouffer_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_MI.mat']), 'fp_combined_fisher_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_POS.mat']), 'fp_combined_fisher_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_z_rayleigh.mat']), 'fp_combined_fisher_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_U2watson.mat']), 'fp_combined_fisher_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_rms_logregress.mat']), 'fp_combined_fisher_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_MI.mat']), 'fp_combined_edington_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_POS.mat']), 'fp_combined_edington_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_z_rayleigh.mat']), 'fp_combined_edington_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_U2watson.mat']), 'fp_combined_edington_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_rms_logregress.mat']), 'fp_combined_edington_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_MI.mat']), 'fp_combined_tippett_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_POS.mat']), 'fp_combined_tippett_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_z_rayleigh.mat']), 'fp_combined_tippett_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_U2watson.mat']), 'fp_combined_tippett_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_rms_logregress.mat']), 'fp_combined_tippett_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_MI.mat']), 'fp_combined_friston_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_POS.mat']), 'fp_combined_friston_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_z_rayleigh.mat']), 'fp_combined_friston_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_U2watson.mat']), 'fp_combined_friston_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_rms_logregress.mat']), 'fp_combined_friston_rms_logregress');

% Sensitivity

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_MI.mat']), 'sensitivity_empirical_vs_chance_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_POS.mat']), 'sensitivity_empirical_vs_chance_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_z_rayleigh.mat']), 'sensitivity_empirical_vs_chance_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_U2watson.mat']), 'sensitivity_empirical_vs_chance_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_rms_logregress.mat']), 'sensitivity_empirical_vs_chance_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_MI.mat']), 'sensitivity_surravg_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_POS.mat']), 'sensitivity_surravg_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_z_rayleigh.mat']), 'sensitivity_surravg_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_U2watson.mat']), 'sensitivity_surravg_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_rms_logregress.mat']), 'sensitivity_surravg_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_MI.mat']), 'sensitivity_combined_stouffer_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_POS.mat']), 'sensitivity_combined_stouffer_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_z_rayleigh.mat']), 'sensitivity_combined_stouffer_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_U2watson.mat']), 'sensitivity_combined_stouffer_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_rms_logregress.mat']), 'sensitivity_combined_stouffer_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_MI.mat']), 'sensitivity_combined_fisher_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_POS.mat']), 'sensitivity_combined_fisher_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_z_rayleigh.mat']), 'sensitivity_combined_fisher_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_U2watson.mat']), 'sensitivity_combined_fisher_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_rms_logregress.mat']), 'sensitivity_combined_fisher_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_MI.mat']), 'sensitivity_combined_edington_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_POS.mat']), 'sensitivity_combined_edington_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_z_rayleigh.mat']), 'sensitivity_combined_edington_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_U2watson.mat']), 'sensitivity_combined_edington_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_rms_logregress.mat']), 'sensitivity_combined_edington_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_MI.mat']), 'sensitivity_combined_tippett_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_POS.mat']), 'sensitivity_combined_tippett_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_z_rayleigh.mat']), 'sensitivity_combined_tippett_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_U2watson.mat']), 'sensitivity_combined_tippett_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_rms_logregress.mat']), 'sensitivity_combined_tippett_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_MI.mat']), 'sensitivity_combined_friston_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_POS.mat']), 'sensitivity_combined_friston_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_z_rayleigh.mat']), 'sensitivity_combined_friston_z_rayleigh');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_U2watson.mat']), 'sensitivity_combined_friston_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_rms_logregress.mat']), 'sensitivity_combined_friston_rms_logregress');

figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w');

% MI
subplot(2,3,1);
stou=plot(sensitivity_combined_stouffer_MI, '-.*', 'LineWidth', 2, 'Color', [0.3 0 0.5]);
hold on;
fish=plot(sensitivity_combined_fisher_MI, ':o', 'LineWidth', 2, 'Color', [0 0.4 0]);
hold on;
edg=plot(sensitivity_combined_edington_MI, '-s', 'LineWidth', 2, 'Color', [0 0.75 1]);
hold on;
sa=plot(sensitivity_surravg_MI, '--v', 'LineWidth', 2, 'Color', [0.5 0 0]);
hold on;
evc=plot(sensitivity_empirical_vs_chance_MI, '-<', 'LineWidth', 2, 'Color', [1 0.27 0]);
xticks(1:8);
xticklabels({'5%','10%','15%','20%','25%','30%','35%','40%'});
yticks(0.2:0.2:1);
ax = gca;
ax.FontSize = 13;
xlim([1 8]);
ylim([0 1]);
xlabel('phase-outcome couplling strength (%)', 'FontSize', 15);
ylabel('% True Positives', 'FontSize', 15);
title('MI', 'FontSize', 17);
% grid on;

% Watson
subplot(2,3,2);
stou=plot(sensitivity_combined_stouffer_U2watson, '-.*', 'LineWidth', 2, 'Color', [0.3 0 0.5]);
hold on;
fish=plot(sensitivity_combined_fisher_U2watson, ':o', 'LineWidth', 2, 'Color', [0 0.4 0]);
hold on;
edg=plot(sensitivity_combined_edington_U2watson, '-s', 'LineWidth', 2, 'Color', [0 0.75 1]);
hold on;
sa=plot(sensitivity_surravg_U2watson, '--v', 'LineWidth', 2, 'Color', [0.5 0 0]);
hold on;
evc=plot(sensitivity_empirical_vs_chance_U2watson, '-<', 'LineWidth', 2, 'Color', [1 0.27 0]);
xticks(1:length(cfg_simulations.strengths_phase_outcome_coupling)-1);
xticklabels({'5%','10%','15%','20%','25%','30%','35%','40%'});
yticks(0.2:0.2:1);
ax = gca;
ax.FontSize = 13;
xlim([1 8]);
ylim([0 1]);
xlabel('phase-outcome couplling strength (%)', 'FontSize', 15);
ylabel('% True Positives', 'FontSize', 15);
title('Watson', 'FontSize', 17);
% grid on;

% Circular logistic regression
subplot(2,3,3);
stou=plot(sensitivity_combined_stouffer_rms_logregress, '-.*', 'LineWidth', 2, 'Color', [0.3 0 0.5]);
hold on;
fish=plot(sensitivity_combined_fisher_rms_logregress, ':o', 'LineWidth', 2, 'Color', [0 0.4 0]);
hold on;
edg=plot(sensitivity_combined_edington_rms_logregress, '-s', 'LineWidth', 2, 'Color', [0 0.75 1]);
hold on;
sa=plot(sensitivity_surravg_rms_logregress, '--v', 'LineWidth', 2, 'Color', [0.5 0 0]);
hold on;
evc=plot(sensitivity_empirical_vs_chance_rms_logregress, '-<', 'LineWidth', 2, 'Color', [1 0.27 0]);
xticks(1:length(cfg_simulations.strengths_phase_outcome_coupling)-1);
xticklabels({'5%','10%','15%','20%','25%','30%','35%','40%'});
yticks(0.2:0.2:1);
lgd=legend([evc sa stou fish edg], 't-test empirical vs. chance', 'surrogate average', 'combined Stouffer', 'combined Fisher', 'combined edington', 'Location','northeast');
lgd.FontSize = 15;
ax = gca;
ax.FontSize = 13;
xlim([1 8]);
ylim([0 1]);
xlabel('phase-outcome couplling strength (%)', 'FontSize', 15);
ylabel('% True Positives', 'FontSize', 15);
title('Circ. log. regr.', 'FontSize', 17);
% grid on;

% POS
subplot(2,3,4);
stou=plot(sensitivity_combined_stouffer_POS, '-.*', 'LineWidth', 2, 'Color', [0.3 0 0.5]);
hold on;
fish=plot(sensitivity_combined_fisher_POS, ':o', 'LineWidth', 2, 'Color', [0 0.4 0]);
hold on;
edg=plot(sensitivity_combined_edington_POS, '-s', 'LineWidth', 2, 'Color', [0 0.75 1]);
hold on;
sa=plot(sensitivity_surravg_POS, '--v', 'LineWidth', 2, 'Color', [0.5 0 0]);
hold on;
evc=plot(sensitivity_empirical_vs_chance_POS, '-<', 'LineWidth', 2, 'Color', [1 0.27 0]);
xticks(1:length(cfg_simulations.strengths_phase_outcome_coupling)-1);
xticklabels({'5%','10%','15%','20%','25%','30%','35%','40%'});
yticks(0.2:0.2:1);
ax = gca;
ax.FontSize = 13;
xlim([1 8]);
ylim([0 1]);
xlabel('phase-outcome couplling strength (%)', 'FontSize', 15);
ylabel('% True Positives', 'FontSize', 15);
title('POS', 'FontSize', 17);
% grid on;

% Rayleigh
subplot(2,3,5);
stou=plot(sensitivity_combined_stouffer_z_rayleigh, '-.*', 'LineWidth', 2, 'Color', [0.3 0 0.5]);
hold on;
fish=plot(sensitivity_combined_fisher_z_rayleigh, ':o', 'LineWidth', 2, 'Color', [0 0.4 0]);
hold on;
edg=plot(sensitivity_combined_edington_z_rayleigh, '-s', 'LineWidth', 2, 'Color', [0 0.75 1]);
hold on;
sa=plot(sensitivity_surravg_z_rayleigh, '--v', 'LineWidth', 2, 'Color', [0.5 0 0]);
hold on;
evc=plot(sensitivity_empirical_vs_chance_z_rayleigh, '-<', 'LineWidth', 2, 'Color', [1 0.27 0]);
xticks(1:length(cfg_simulations.strengths_phase_outcome_coupling)-1);
xticklabels({'5%','10%','15%','20%','25%','30%','35%','40%'});
yticks(0.2:0.2:1);
ax = gca;
ax.FontSize = 13;
xlim([1 8]);
ylim([0 1]);
xlabel('phase-outcome couplling strength (%)', 'FontSize', 15);
ylabel('% True Positives', 'FontSize', 15);
title('Rayleigh', 'FontSize', 17);
% grid on;

% show FP rate for each method
figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w');
subplot(2,3,1);
evc=bar(1, fp_empirical_vs_chance_MI); set(evc,'FaceColor', [1 0.27 0]); hold on;
sa=bar(2, fp_surravg_MI); set(sa,'FaceColor',[0.5 0 0]); hold on;
stou=bar(3, fp_combined_stouffer_MI); set(stou,'FaceColor',[0.3 0 0.5]); hold on;
fish=bar(4, fp_combined_fisher_MI); set(fish,'FaceColor',[0 0.4 0]); hold on;
edg=bar(5, fp_combined_edington_MI); set(edg,'FaceColor',[0 0.75 1]); hold on;
hline(0.05, '--k');
% ylim([0 0.1]);
xticks(1:6);
xticklabels({'t-test empirical vs. chance','surrogate average','combined Stouffer','combined Fisher','combined Edington'});
xtickangle(45);
ylabel('False Positive rate (%)', 'FontSize', 15);
title('MI', 'FontSize', 15);
ax = gca;
ax.FontSize = 12;
subplot(2,3,2);
evc=bar(1, fp_empirical_vs_chance_U2watson); set(evc,'FaceColor', [1 0.27 0]); hold on;
sa=bar(2, fp_surravg_U2watson); set(sa,'FaceColor',[0.5 0 0]); hold on;
stou=bar(3, fp_combined_stouffer_U2watson); set(stou,'FaceColor',[0.3 0 0.5]); hold on;
fish=bar(4, fp_combined_fisher_U2watson); set(fish,'FaceColor',[0 0.4 0]); hold on;
edg=bar(5, fp_combined_edington_U2watson); set(edg,'FaceColor',[0 0.75 1]); hold on;
hline(0.05, '--k');
% ylim([0 0.1]);
xticks(1:6);
xticklabels({'t-test empirical vs. chance','surrogate average','combined Stouffer','combined Fisher','combined Edington'});
xtickangle(45);
ylabel('False Positive rate (%)', 'FontSize', 15);
title('POS', 'FontSize', 15);
ax = gca;
ax.FontSize = 12;
subplot(2,3,3);
evc=bar(1, fp_empirical_vs_chance_rms_logregress); set(evc,'FaceColor', [1 0.27 0]); hold on;
sa=bar(2, fp_surravg_rms_logregress); set(sa,'FaceColor',[0.5 0 0]); hold on;
stou=bar(3, fp_combined_stouffer_rms_logregress); set(stou,'FaceColor',[0.3 0 0.5]); hold on;
fish=bar(4, fp_combined_fisher_rms_logregress); set(fish,'FaceColor',[0 0.4 0]); hold on;
edg=bar(5, fp_combined_edington_rms_logregress); set(edg,'FaceColor',[0 0.75 1]); hold on;
hline(0.05, '--k');
% ylim([0 0.1]);
xticks(1:6);
xticklabels({'t-test empirical vs. chance','surrogate average','combined Stouffer','combined Fisher','combined Edington'});
xtickangle(45);
ylabel('False Positive rate (%)', 'FontSize', 15);
title('Rayleigh', 'FontSize', 15);
ax = gca;
ax.FontSize = 12;
subplot(2,3,4);
evc=bar(1, fp_empirical_vs_chance_POS); set(evc,'FaceColor', [1 0.27 0]); hold on;
sa=bar(2, fp_surravg_POS); set(sa,'FaceColor',[0.5 0 0]); hold on;
stou=bar(3, fp_combined_stouffer_POS); set(stou,'FaceColor',[0.3 0 0.5]); hold on;
fish=bar(4, fp_combined_fisher_POS); set(fish,'FaceColor',[0 0.4 0]); hold on;
edg=bar(5, fp_combined_edington_POS); set(edg,'FaceColor',[0 0.75 1]); hold on;
hline(0.05, '--k');
% ylim([0 0.1]);
xticks(1:6);
xticklabels({'t-test empirical vs. chance','surrogate average','combined Stouffer','combined Fisher','combined Edington'});
xtickangle(45);
ylabel('False Positive rate (%)', 'FontSize', 15);
title('Watson', 'FontSize', 15);
ax = gca;
ax.FontSize = 12;
subplot(2,3,5);
evc=bar(1, fp_empirical_vs_chance_z_rayleigh); set(evc,'FaceColor', [1 0.27 0]); hold on;
sa=bar(2, fp_surravg_z_rayleigh); set(sa,'FaceColor',[0.5 0 0]); hold on;
stou=bar(3, fp_combined_stouffer_z_rayleigh); set(stou,'FaceColor',[0.3 0 0.5]); hold on;
fish=bar(4, fp_combined_fisher_z_rayleigh); set(fish,'FaceColor',[0 0.4 0]); hold on;
edg=bar(5, fp_combined_edington_z_rayleigh); set(edg,'FaceColor',[0 0.75 1]); hold on;
hline(0.05, '--k');
% ylim([0 0.1]);
xticks(1:6);
xticklabels({'t-test empirical vs. chance','surrogate average','combined Stouffer','combined Fisher','combined Edington'});
xtickangle(45);
ylabel('False Positive rate (%)', 'FontSize', 15);
title('Log. regr.', 'FontSize', 15);
ax = gca;
ax.FontSize = 12;

end

