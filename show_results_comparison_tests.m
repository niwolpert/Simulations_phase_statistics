function show_results_comparison_tests(cfg_simulations, root_dir, method_estimate_significance)
%SHOW_RESULTS_COMPARISON_TESTS
%   Shows the results to compare the sensitivity and False positive rate 
%   of the five circular statistical tests.
%   
%   INPUTS
%   - cfg_simulations:                 Configuration structure with
%                                      simulation parameters
%   - root_dir:                        Root directory where results will be
%                                      saved
%   - method_estimate_significance:    Method to estimate significance on
%                                      the group level (string).
%                                      Options:
%                                      - 'empirical_vs_chance'
%                                      - 'surravg'
%                                      - 'combined_stouffer'
%                                      - 'combined_fisher'
%                                      - 'combined_edington'
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
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_U2watson.mat']), 'fp_empirical_vs_chance_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_rms_logregress.mat']), 'fp_empirical_vs_chance_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_MI.mat']), 'fp_surravg_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_POS.mat']), 'fp_surravg_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_U2watson.mat']), 'fp_surravg_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_rms_logregress.mat']), 'fp_surravg_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_MI.mat']), 'fp_combined_stouffer_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_POS.mat']), 'fp_combined_stouffer_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_U2watson.mat']), 'fp_combined_stouffer_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_rms_logregress.mat']), 'fp_combined_stouffer_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_MI.mat']), 'fp_combined_fisher_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_POS.mat']), 'fp_combined_fisher_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_U2watson.mat']), 'fp_combined_fisher_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_rms_logregress.mat']), 'fp_combined_fisher_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_MI.mat']), 'fp_combined_edington_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_POS.mat']), 'fp_combined_edington_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_U2watson.mat']), 'fp_combined_edington_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_rms_logregress.mat']), 'fp_combined_edington_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_MI.mat']), 'fp_combined_tippett_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_POS.mat']), 'fp_combined_tippett_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_U2watson.mat']), 'fp_combined_tippett_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_rms_logregress.mat']), 'fp_combined_tippett_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_MI.mat']), 'fp_combined_friston_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_POS.mat']), 'fp_combined_friston_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_U2watson.mat']), 'fp_combined_friston_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_rms_logregress.mat']), 'fp_combined_friston_rms_logregress');

% Sensitivity

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_MI.mat']), 'sensitivity_empirical_vs_chance_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_POS.mat']), 'sensitivity_empirical_vs_chance_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_U2watson.mat']), 'sensitivity_empirical_vs_chance_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_rms_logregress.mat']), 'sensitivity_empirical_vs_chance_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_MI.mat']), 'sensitivity_surravg_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_POS.mat']), 'sensitivity_surravg_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_U2watson.mat']), 'sensitivity_surravg_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_rms_logregress.mat']), 'sensitivity_surravg_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_MI.mat']), 'sensitivity_combined_stouffer_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_POS.mat']), 'sensitivity_combined_stouffer_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_U2watson.mat']), 'sensitivity_combined_stouffer_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_rms_logregress.mat']), 'sensitivity_combined_stouffer_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_MI.mat']), 'sensitivity_combined_fisher_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_POS.mat']), 'sensitivity_combined_fisher_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_U2watson.mat']), 'sensitivity_combined_fisher_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_rms_logregress.mat']), 'sensitivity_combined_fisher_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_MI.mat']), 'sensitivity_combined_edington_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_POS.mat']), 'sensitivity_combined_edington_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_U2watson.mat']), 'sensitivity_combined_edington_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_rms_logregress.mat']), 'sensitivity_combined_edington_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_MI.mat']), 'sensitivity_combined_tippett_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_POS.mat']), 'sensitivity_combined_tippett_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_U2watson.mat']), 'sensitivity_combined_tippett_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_rms_logregress.mat']), 'sensitivity_combined_tippett_rms_logregress');

load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_MI.mat']), 'sensitivity_combined_friston_MI');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_POS.mat']), 'sensitivity_combined_friston_POS');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_U2watson.mat']), 'sensitivity_combined_friston_U2watson');
load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_rms_logregress.mat']), 'sensitivity_combined_friston_rms_logregress');


if strcmp(method_estimate_significance, 'empirical_vs_chance')
    fp_rms_logregress = fp_empirical_vs_chance_rms_logregress;
    fp_POS = fp_empirical_vs_chance_POS;
    fp_U2watson = fp_empirical_vs_chance_U2watson;
    fp_MI = fp_empirical_vs_chance_MI;
    
    sensitivity_rms_logregress = sensitivity_empirical_vs_chance_rms_logregress;
    sensitivity_POS = sensitivity_empirical_vs_chance_POS;
    sensitivity_U2watson = sensitivity_empirical_vs_chance_U2watson;
    sensitivity_MI = sensitivity_empirical_vs_chance_MI;

elseif strcmp(method_estimate_significance, 'surravg')
    fp_rms_logregress = fp_surravg_rms_logregress;
    fp_POS = fp_surravg_POS;
    fp_U2watson = fp_surravg_U2watson;
    fp_MI = fp_surravg_MI;
    
    sensitivity_rms_logregress = sensitivity_surravg_rms_logregress;
    sensitivity_POS = sensitivity_surravg_POS;
    sensitivity_U2watson = sensitivity_surravg_U2watson;
    sensitivity_MI = sensitivity_surravg_MI;
    
elseif strcmp(method_estimate_significance, 'combined_stouffer')
    fp_rms_logregress = fp_combined_stouffer_rms_logregress;
    fp_POS = fp_combined_stouffer_POS;
    fp_U2watson = fp_combined_stouffer_U2watson;
    fp_MI = fp_combined_stouffer_MI;
    
    sensitivity_rms_logregress = sensitivity_combined_stouffer_rms_logregress;
    sensitivity_POS = sensitivity_combined_stouffer_POS;
    sensitivity_U2watson = sensitivity_combined_stouffer_U2watson;
    sensitivity_MI = sensitivity_combined_stouffer_MI;
    
elseif strcmp(method_estimate_significance, 'combined_fisher')
    fp_rms_logregress = fp_combined_fisher_rms_logregress;
    fp_POS = fp_combined_fisher_POS;
    fp_U2watson = fp_combined_fisher_U2watson;
    fp_MI = fp_combined_fisher_MI;
    
    sensitivity_rms_logregress = sensitivity_combined_fisher_rms_logregress;
    sensitivity_POS = sensitivity_combined_fisher_POS;
    sensitivity_U2watson = sensitivity_combined_fisher_U2watson;
    sensitivity_MI = sensitivity_combined_fisher_MI;
    
elseif strcmp(method_estimate_significance, 'combined_edington')
    fp_rms_logregress = fp_combined_edington_rms_logregress;
    fp_POS = fp_combined_edington_POS;
    fp_U2watson = fp_combined_edington_U2watson;
    fp_MI = fp_combined_edington_MI;
    
    sensitivity_rms_logregress = sensitivity_combined_edington_rms_logregress;
    sensitivity_POS = sensitivity_combined_edington_POS;
    sensitivity_U2watson = sensitivity_combined_edington_U2watson;
    sensitivity_MI = sensitivity_combined_edington_MI;
end

close all;

colors = {'r',[1 0.65 0],[0 0.5 0],[0.12 0.56 1]};

figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','w');
subplot(1,3,1:2);
l=plot(sensitivity_rms_logregress(2:end), 'Color', colors{1}, 'Linewidth', 3); hold on;
p=plot(sensitivity_POS(2:end), 'Color', colors{2}, 'Linewidth', 3); hold on;
u=plot(sensitivity_U2watson(2:end), 'Color', colors{3}, 'Linewidth', 3); hold on;
m=plot(sensitivity_MI(2:end),'Color',  colors{4}, 'Linewidth', 3); hold on;
xticks(1:8);
xticklabels({'5','10','15','20','25','30','35','40'});
lgd=legend([m p u l], 'MI', 'POS', 'Watson', 'circ. log. regr.', 'Location','northeast');
lgd.FontSize = 15;
ax = gca;
ax.FontSize = 16;
xlim([1 8]);
ylim([0 1]);
xlabel('phase-outcome coupling strength', 'FontSize', 20);
ylabel('% True Positives', 'FontSize', 20);


if strcmp(method_estimate_significance, 'empirical_vs_chance')

    title(['Sensitivity, t-test empirical vs. chance, ' cfg_simulations.oscillation ', coupling mode ' num2str(cfg_simulations.coupling_mode) ':1, ' num2str(cfg_simulations.nTrials_poc) ' trials'], 'FontSize', 20, 'Interpreter', 'none');
    
elseif strcmp(method_estimate_significance, 'surravg')

    title(['Sensitivity, surrogate average, ' cfg_simulations.oscillation ', coupling mode ' num2str(cfg_simulations.coupling_mode) ':1, ' num2str(cfg_simulations.nTrials_poc) ' trials'], 'FontSize', 20, 'Interpreter', 'none');
    
elseif strcmp(method_estimate_significance, 'combined_stouffer')

    title(['Sensitivity, combined stouffer, ' cfg_simulations.oscillation ', coupling mode ' num2str(cfg_simulations.coupling_mode) ':1, ' num2str(cfg_simulations.nTrials_poc) ' trials'], 'FontSize', 20, 'Interpreter', 'none');
    
elseif strcmp(method_estimate_significance, 'combined_fisher')

    title(['Sensitivity, combined fisher, ' cfg_simulations.oscillation ', coupling mode ' num2str(cfg_simulations.coupling_mode) ':1, ' num2str(cfg_simulations.nTrials_poc) ' trials'], 'FontSize', 20, 'Interpreter', 'none');
    
elseif strcmp(method_estimate_significance, 'combined_edington')
    
    title(['Sensitivity, combined edington, ' cfg_simulations.oscillation ', coupling mode ' num2str(cfg_simulations.coupling_mode) ':1, ' num2str(cfg_simulations.nTrials_poc) ' trials'], 'FontSize', 20, 'Interpreter', 'none');
    
end

grid on;

subplot(1,3,3);
l=bar(1, fp_rms_logregress); set(l,'FaceColor', colors{1}); hold on;
p=bar(2, fp_POS); set(p,'FaceColor',colors{2}); hold on;
u=bar(3, fp_U2watson); set(u,'FaceColor', colors{3}); hold on;
m=bar(4, fp_MI); set(m,'FaceColor', colors{4}); hold on;
hline(0.05, '--k');
xlim([0.4 4.6]);
ylim([0 0.06]);
xticks(1:4);
xticklabels({'Cir. log. regr.','POS','Watson','MI'});
xtickangle(45);
ylabel('% False Positives', 'FontSize', 20);
title(['FP-rate, ' method_estimate_significance], 'FontSize', 20, 'interpreter', 'none');
ax = gca;
ax.FontSize = 14;
% grid on;

end

