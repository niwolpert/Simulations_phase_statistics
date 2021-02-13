function show_results_amplitude(cfg_simulations, root_dir)
%SHOW_RESULTS_AMOPLITUDE
%   Shows the results of simulations on amplitude.
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

load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep 'amplitudes_sensitivity_coupling_strength_ttest_MI.mat']), 'sensitivity_coupling_strength_ttest_MI');
load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep 'amplitudes_sensitivity_coupling_strength_ttest_POS.mat']), 'sensitivity_coupling_strength_ttest_POS');
load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep 'amplitudes_sensitivity_coupling_strength_ttest_U2watson.mat']), 'sensitivity_coupling_strength_ttest_U2watson');
load(strcat([root_dir 'Amplitude' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep 'amplitudes_sensitivity_coupling_strength_ttest_rms_logregress.mat']), 'sensitivity_coupling_strength_ttest_rms_logregress');

close all;

colors = {'r',[1 0.65 0],[0 0.5 0],[0.12 0.56 1]};

figure('units','normalized','outerposition',[0.1105 0.2229 0.5590 0.6083]); set(gcf,'color','w');
l=plot(cfg_simulations.steps_amplitudes, sensitivity_coupling_strength_ttest_rms_logregress*100, 'Color', colors{1}, 'Linewidth', 3); hold on;
p=plot(cfg_simulations.steps_amplitudes, sensitivity_coupling_strength_ttest_POS*100, 'Color', colors{2}, 'Linewidth', 3); hold on;
u=plot(cfg_simulations.steps_amplitudes, sensitivity_coupling_strength_ttest_U2watson*100, 'Color', colors{3}, 'Linewidth', 3); hold on;
m=plot(cfg_simulations.steps_amplitudes, sensitivity_coupling_strength_ttest_MI*100,'Color',  colors{4}, 'Linewidth', 3); hold on;
hline(5, '--k');
lgd=legend([m p u l], 'MI', 'POS', 'Watson', 'circ. log. regr.', 'Location','northeast');
lgd.FontSize = 15;
ax = gca;
ax.FontSize = 16;
xlabel('Amplitude', 'FontSize', 20);
ylabel('% True Positives', 'FontSize', 20);
grid on;

end
