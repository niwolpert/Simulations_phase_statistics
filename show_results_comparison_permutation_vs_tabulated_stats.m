function show_results_comparison_permutation_vs_tabulated_stats(cfg_simulations, root_dir)
%SHOW_RESULTS_COMPARISON_PERMUTATION_VS_TABULATED_STATS
%   Shows the results to compares the sensitivity and False positive rate 
%   for the Rayleigh test, Watson test and circular logistic regression
%   using permutations vs. tabulated statistics.
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

% merge p-values from single-subject level
rayleigh_pvalues_direct_all = [];
watson_pvalues_direct_all = [];
logregress_pvalues_direct_all = [];
rayleigh_pvalues_permutation_all = [];
watson_pvalues_permutation_all = [];
logregress_pvalues_permutation_all = [];

dispstat('','init');
for iexp=1:1000
    
    dispstat(sprintf('Loading stats from experiment %d%',iexp));
    
    if isfile(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_experiment' num2str(iexp) '.mat']))
        
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_empirical_experiment' num2str(iexp) '.mat']), 'z_rayleigh_empirical');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_experiment' num2str(iexp) '.mat']), 'U2watson_empirical');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical');
        
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_surr_distr_experiment' num2str(iexp) '.mat']), 'z_rayleigh_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr');
        
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues_direct' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rayleigh_pvalues_direct_experiment' num2str(iexp) '.mat']), 'rayleigh_pvalues_direct');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues_direct' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_watson_pvalues_direct_experiment' num2str(iexp) '.mat']), 'watson_pvalues_direct');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues_direct' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_logregress_pvalues_direct_experiment' num2str(iexp) '.mat']), 'logregress_pvalues_direct');
        
        % note direct p-values
        rayleigh_pvalues_direct_all = [rayleigh_pvalues_direct_all; rayleigh_pvalues_direct];
        watson_pvalues_direct_all = [watson_pvalues_direct_all; watson_pvalues_direct];
        logregress_pvalues_direct_all = [logregress_pvalues_direct_all; logregress_pvalues_direct];
        
        % compute p-value as percentage of surrogates higher than empirical
        rayleigh_pvalues_permutation = nan(30, length(cfg_simulations.strengths_phase_outcome_coupling));
        watson_pvalues_permutation = nan(30, length(cfg_simulations.strengths_phase_outcome_coupling));
        logregress_pvalues_permutation = nan(30, length(cfg_simulations.strengths_phase_outcome_coupling));
        for isubject=1:30
            
            for ipoc=1:length(cfg_simulations.strengths_phase_outcome_coupling)
                
                rayleigh_pvalues_permutation(isubject, ipoc) = length(find(z_rayleigh_surr_distr{isubject, ipoc}>z_rayleigh_empirical(isubject, ipoc)))/length(find(~isnan(z_rayleigh_surr_distr{isubject, ipoc})));
                watson_pvalues_permutation(isubject, ipoc) = length(find(U2watson_surr_distr{isubject, ipoc}>U2watson_empirical(isubject, ipoc)))/length(find(~isnan(U2watson_surr_distr{isubject, ipoc})));
                logregress_pvalues_permutation(isubject, ipoc) = length(find(rms_logregress_surr_distr{isubject, ipoc}>rms_logregress_empirical(isubject, ipoc)))/length(find(~isnan(rms_logregress_surr_distr{isubject, ipoc})));
                
            end
        end
        rayleigh_pvalues_permutation_all = [rayleigh_pvalues_permutation_all; rayleigh_pvalues_permutation];
        watson_pvalues_permutation_all = [watson_pvalues_permutation_all; watson_pvalues_permutation];
        logregress_pvalues_permutation_all = [logregress_pvalues_permutation_all; logregress_pvalues_permutation];
    end
end

% compare sensitivity for direct vs. permutation
rayleigh_sensitivity_direct = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
watson_sensitivity_direct = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
logregress_sensitivity_direct = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
rayleigh_sensitivity_permutation = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
watson_sensitivity_permutation = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
logregress_sensitivity_permutation = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);

thresh_signif = 0.05;
for ipoc=2:length(cfg_simulations.strengths_phase_outcome_coupling)
    
    rayleigh_sensitivity_direct(ipoc-1) = length(find(rayleigh_pvalues_direct_all(:, ipoc)<thresh_signif))/length(find(~isnan(rayleigh_pvalues_direct_all(:, ipoc))));
    watson_sensitivity_direct(ipoc-1) = length(find(watson_pvalues_direct_all(:, ipoc)<thresh_signif))/length(find(~isnan(watson_pvalues_direct_all(:, ipoc))));
    logregress_sensitivity_direct(ipoc-1) = length(find(logregress_pvalues_direct_all(:, ipoc)<thresh_signif))/length(find(~isnan(logregress_pvalues_direct_all(:, ipoc))));
    rayleigh_sensitivity_permutation(ipoc-1) = length(find(rayleigh_pvalues_permutation_all(:, ipoc)<thresh_signif))/length(find(~isnan(rayleigh_pvalues_permutation_all(:, ipoc))));
    watson_sensitivity_permutation(ipoc-1) = length(find(watson_pvalues_permutation_all(:, ipoc)<thresh_signif))/length(find(~isnan(watson_pvalues_permutation_all(:, ipoc))));
    logregress_sensitivity_permutation(ipoc-1) = length(find(logregress_pvalues_permutation_all(:, ipoc)<thresh_signif))/length(find(~isnan(logregress_pvalues_permutation_all(:, ipoc))));
    
end
rayleigh_fp_direct = length(find(rayleigh_pvalues_direct_all(:, 1)<thresh_signif))/length(find(~isnan(rayleigh_pvalues_direct_all(:, 1))));
watson_fp_direct = length(find(watson_pvalues_direct_all(:, 1)<thresh_signif))/length(find(~isnan(watson_pvalues_direct_all(:, 1))));
logregress_fp_direct = length(find(logregress_pvalues_direct_all(:, 1)<thresh_signif))/length(find(~isnan(logregress_pvalues_direct_all(:, 1))));

rayleigh_fp_permutation = length(find(rayleigh_pvalues_permutation_all(:, 1)<thresh_signif))/length(find(~isnan(rayleigh_pvalues_permutation_all(:, 1))));
watson_fp_permutation = length(find(watson_pvalues_permutation_all(:, 1)<thresh_signif))/length(find(~isnan(watson_pvalues_permutation_all(:, 1))));
logregress_fp_permutation = length(find(logregress_pvalues_permutation_all(:, 1)<thresh_signif))/length(find(~isnan(logregress_pvalues_permutation_all(:, 1))));

close all;

figure('units','normalized','outerposition',[0.1055 0.2236 0.6965 0.5750]); set(gcf,'color','w');

subplot(1,3,1:2);
p=plot(rayleigh_sensitivity_permutation, 'Color', 'k', 'Linewidth', 3); hold on;
d=plot(rayleigh_sensitivity_direct, 'Color', [0.7 0.7 0.7], 'Linewidth', 3); hold on;
xticks(1:8);
xticklabels({'5%','10%','15%','20%','25%','30%','35%','40%'});
lgd=legend([p d], 'permutations', 'direct output', 'Location','northeast');
lgd.FontSize = 15;
ax = gca;
ax.FontSize = 16;
xlim([1 8]);
ylim([0 0.7]);
xlabel('phase-locking strength', 'FontSize', 20);
ylabel('Sensitivity', 'FontSize', 20);
title(['Sensitivity Rayleigh, coupling mode ' num2str(cfg_simulations.coupling_mode) ':1'], 'FontSize', 25, 'Interpreter', 'none');
% grid on;

subplot(1,3,3);
p=bar(1, rayleigh_fp_permutation); set(p,'FaceColor', 'k'); hold on;
d=bar(2, rayleigh_fp_direct); set(d,'FaceColor', [0.7 0.7 0.7]); hold on;
hline(0.05, '--k');
xlim([0.4 2.6]);
ylim([0 0.1]);
xticks(1:2);
xticklabels({'permutations', 'direct output'});
xtickangle(45);
yticks(0:0.01:0.1);
ylabel('FP rate (%)', 'FontSize', 20);
title('FP-rate', 'FontSize', 25, 'interpreter', 'none');
ax = gca;
ax.FontSize = 16;
% grid on;

figure('units','normalized','outerposition',[0.1055 0.2236 0.6965 0.5750]); set(gcf,'color','w');

subplot(1,3,1:2);
p=plot(watson_sensitivity_permutation, 'Color', 'k', 'Linewidth', 3); hold on;
d=plot(watson_sensitivity_direct, 'Color', [0.7 0.7 0.7], 'Linewidth', 3); hold on;
xticks(1:8);
xticklabels({'5%','10%','15%','20%','25%','30%','35%','40%'});
lgd=legend([p d], 'permutations', 'direct output', 'Location','northeast');
lgd.FontSize = 15;
ax = gca;
ax.FontSize = 16;
xlim([1 8]);
ylim([0 0.7]);
xlabel('phase-locking strength', 'FontSize', 20);
ylabel('Sensitivity', 'FontSize', 25);
title(['Sensitivity Watson, coupling mode ' num2str(cfg_simulations.coupling_mode) ':1'], 'FontSize', 20, 'Interpreter', 'none');
% grid on;

subplot(1,3,3);
p=bar(1, watson_fp_permutation); set(p,'FaceColor', 'k'); hold on;
d=bar(2, watson_fp_direct); set(d,'FaceColor', [0.7 0.7 0.7]); hold on;
hline(0.05, '--k');
xlim([0.4 2.6]);
ylim([0 0.1]);
xticks(1:2);
xticklabels({'permutations', 'direct output'});
xtickangle(45);
yticks(0:0.01:0.1);
ylabel('FP rate (%)', 'FontSize', 20);
title('FP-rate', 'FontSize', 25, 'interpreter', 'none');
ax = gca;
ax.FontSize = 16;
% grid on;

figure('units','normalized','outerposition',[0.1055 0.2236 0.6965 0.5750]); set(gcf,'color','w');

subplot(1,3,1:2);
p=plot(logregress_sensitivity_permutation, 'Color', 'k', 'Linewidth', 3); hold on;
d=plot(logregress_sensitivity_direct, 'Color', [0.7 0.7 0.7], 'Linewidth', 3); hold on;
xticks(1:8);
xticklabels({'5%','10%','15%','20%','25%','30%','35%','40%'});
lgd=legend([p d], 'permutations', 'direct output', 'Location','northeast');
lgd.FontSize = 15;
ax = gca;
ax.FontSize = 16;
xlim([1 8]);
ylim([0 0.7]);
xlabel('phase-locking strength', 'FontSize', 20);
ylabel('Sensitivity', 'FontSize', 25);
title(['Sensitivity Logistic regression, coupling mode ' num2str(cfg_simulations.coupling_mode) ':1'], 'FontSize', 20, 'Interpreter', 'none');
% grid on;

subplot(1,3,3);
p=bar(1, logregress_fp_permutation); set(p,'FaceColor', 'k'); hold on;
d=bar(2, logregress_fp_direct); set(d,'FaceColor', [0.7 0.7 0.7]); hold on;
hline(0.05, '--k');
xlim([0.4 2.6]);
ylim([0 0.1]);
xticks(1:2);
xticklabels({'permutations', 'direct output'});
xtickangle(45);
yticks(0:0.01:0.1);
ylabel('FP rate (%)', 'FontSize', 20);
title('FP-rate', 'FontSize', 25, 'interpreter', 'none');
ax = gca;
ax.FontSize = 16;
% grid on;

end