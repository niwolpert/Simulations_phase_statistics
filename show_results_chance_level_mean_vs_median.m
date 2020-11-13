function show_results_chance_level_mean_vs_median(cfg_simulations, root_dir)
%SHOW_RESULTS_CHANCE_LEVEL_MEAN_VS_MEDIAN
%   Compares sensitivity and False Positive rate when defining chance level
%   as the mean vs. the median of surrogate distributions
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

% Initialize matrices for noting p-values per virtual experiment and level
% of phase-outcome coupling strength
pvalues_coupling_strength_mean_MI = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_coupling_strength_mean_POS = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_coupling_strength_mean_z_rayleigh = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_coupling_strength_mean_U2watson = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_coupling_strength_mean_rms_logregress = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));

pvalues_coupling_strength_median_MI = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_coupling_strength_median_POS = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_coupling_strength_median_z_rayleigh = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_coupling_strength_median_U2watson = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_coupling_strength_median_rms_logregress = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));

nsubjects = 30;

dispstat('','init');
for iexp=1:nexperiments
    
    dispstat(sprintf('Experiment %d%',iexp));
    
    if isfile(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_experiment' num2str(iexp) '.mat']))
        
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_experiment' num2str(iexp) '.mat']), 'MIs_empirical');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_experiment' num2str(iexp) '.mat']), 'POS_empirical');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_empirical_experiment' num2str(iexp) '.mat']), 'z_rayleigh_empirical');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_experiment' num2str(iexp) '.mat']), 'U2watson_empirical');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical');

        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_experiment' num2str(iexp) '.mat']), 'POS_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_surr_distr_experiment' num2str(iexp) '.mat']), 'z_rayleigh_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr');
        
        for ipoc=1:length(cfg_simulations.strengths_phase_outcome_coupling)
            
            % compute chance levels using mean or median of surrogate
            % distributions
            chance_levels_mean_MI = nan(nsubjects, 1);
            chance_levels_mean_POS = nan(nsubjects, 1);
            chance_levels_mean_rayleigh = nan(nsubjects, 1);
            chance_levels_mean_watson = nan(nsubjects, 1);
            chance_levels_mean_logregress = nan(nsubjects, 1);
            
            chance_levels_median_MI = nan(nsubjects, 1);
            chance_levels_median_POS = nan(nsubjects, 1);
            chance_levels_median_rayleigh = nan(nsubjects, 1);
            chance_levels_median_watson = nan(nsubjects, 1);
            chance_levels_median_logregress = nan(nsubjects, 1);
            for isubject=1:nsubjects
                
                chance_levels_mean_MI(isubject) = nanmean(MIs_surr_distr{isubject, ipoc});
                chance_levels_mean_POS(isubject) = nanmean(POS_surr_distr{isubject, ipoc});
                chance_levels_mean_rayleigh(isubject) = nanmean(z_rayleigh_surr_distr{isubject, ipoc});
                chance_levels_mean_watson(isubject) = nanmean(U2watson_surr_distr{isubject, ipoc});
                chance_levels_mean_logregress(isubject) = nanmean(rms_logregress_surr_distr{isubject, ipoc});
                
                chance_levels_median_MI(isubject) = nanmedian(MIs_surr_distr{isubject, ipoc});
                chance_levels_median_POS(isubject) = nanmedian(POS_surr_distr{isubject, ipoc});
                chance_levels_median_rayleigh(isubject) = nanmedian(z_rayleigh_surr_distr{isubject, ipoc});
                chance_levels_median_watson(isubject) = nanmedian(U2watson_surr_distr{isubject, ipoc});
                chance_levels_median_logregress(isubject) = nanmedian(rms_logregress_surr_distr{isubject, ipoc});
                
            end
            
            % Compute p-values using either mean or median as chance level
            [~,pvalues_coupling_strength_mean_MI(iexp, ipoc)] = ttest(MIs_empirical(:, ipoc), chance_levels_mean_MI,'Tail','right');
            [~,pvalues_coupling_strength_mean_POS(iexp, ipoc)] = ttest(POS_empirical(:, ipoc), chance_levels_mean_POS,'Tail','right');
            [~,pvalues_coupling_strength_mean_z_rayleigh(iexp, ipoc)] = ttest(z_rayleigh_empirical(:, ipoc), chance_levels_mean_rayleigh,'Tail','right');
            [~,pvalues_coupling_strength_mean_U2watson(iexp, ipoc)] = ttest(U2watson_empirical(:, ipoc), chance_levels_mean_watson,'Tail','right');
            [~,pvalues_coupling_strength_mean_rms_logregress(iexp, ipoc)] = ttest(rms_logregress_empirical(:, ipoc), chance_levels_mean_logregress,'Tail','right');
            
            [~,pvalues_coupling_strength_median_MI(iexp, ipoc)] = ttest(MIs_empirical(:, ipoc), chance_levels_median_MI,'Tail','right');
            [~,pvalues_coupling_strength_median_POS(iexp, ipoc)] = ttest(POS_empirical(:, ipoc), chance_levels_median_POS,'Tail','right');
            [~,pvalues_coupling_strength_median_z_rayleigh(iexp, ipoc)] = ttest(z_rayleigh_empirical(:, ipoc), chance_levels_median_rayleigh(isubject),'Tail','right');
            [~,pvalues_coupling_strength_median_U2watson(iexp, ipoc)] = ttest(U2watson_empirical(:, ipoc), chance_levels_median_watson,'Tail','right');
            [~,pvalues_coupling_strength_median_rms_logregress(iexp, ipoc)] = ttest(rms_logregress_empirical(:, ipoc), chance_levels_median_logregress,'Tail','right');
            
        end
    end
end

signific_thresh = 0.05;

% Compute FP-rate

fp_coupling_strength_median_MI = length(find(pvalues_coupling_strength_median_MI(:, 1)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_median_MI(:, 1))));
fp_coupling_strength_median_POS = length(find(pvalues_coupling_strength_median_POS(:, 1)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_median_POS(:, 1))));
fp_coupling_strength_median_z_rayleigh = length(find(pvalues_coupling_strength_median_z_rayleigh(:, 1)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_median_z_rayleigh(:, 1))));
fp_coupling_strength_median_U2watson = length(find(pvalues_coupling_strength_median_U2watson(:, 1)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_median_U2watson(:, 1))));
fp_coupling_strength_median_rms_logregress = length(find(pvalues_coupling_strength_median_rms_logregress(:, 1)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_median_rms_logregress(:, 1))));

fp_coupling_strength_mean_MI = length(find(pvalues_coupling_strength_mean_MI(:, 1)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_mean_MI(:, 1))));
fp_coupling_strength_mean_POS = length(find(pvalues_coupling_strength_mean_POS(:, 1)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_mean_POS(:, 1))));
fp_coupling_strength_mean_z_rayleigh = length(find(pvalues_coupling_strength_mean_z_rayleigh(:, 1)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_mean_z_rayleigh(:, 1))));
fp_coupling_strength_mean_U2watson = length(find(pvalues_coupling_strength_mean_U2watson(:, 1)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_mean_U2watson(:, 1))));
fp_coupling_strength_mean_rms_logregress = length(find(pvalues_coupling_strength_mean_rms_logregress(:, 1)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_mean_rms_logregress(:, 1))));

% Compute sensitivity

sensitivity_coupling_strength_median_MI = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
sensitivity_coupling_strength_median_POS = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
sensitivity_coupling_strength_median_z_rayleigh = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
sensitivity_coupling_strength_median_U2watson = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
sensitivity_coupling_strength_median_rms_logregress = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);

sensitivity_coupling_strength_mean_MI = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
sensitivity_coupling_strength_mean_POS = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
sensitivity_coupling_strength_mean_z_rayleigh = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
sensitivity_coupling_strength_mean_U2watson = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
sensitivity_coupling_strength_mean_rms_logregress = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling)-1);
for ipoc=1:length(cfg_simulations.strengths_phase_outcome_coupling)
    
    sensitivity_coupling_strength_median_MI(ipoc) = length(find(pvalues_coupling_strength_median_MI(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_median_MI(:, ipoc))));
    sensitivity_coupling_strength_median_POS(ipoc) = length(find(pvalues_coupling_strength_median_POS(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_median_POS(:, ipoc))));
    sensitivity_coupling_strength_median_z_rayleigh(ipoc) = length(find(pvalues_coupling_strength_median_z_rayleigh(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_median_z_rayleigh(:, ipoc))));
    sensitivity_coupling_strength_median_U2watson(ipoc) = length(find(pvalues_coupling_strength_median_U2watson(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_median_U2watson(:, ipoc))));
    sensitivity_coupling_strength_median_rms_logregress(ipoc) = length(find(pvalues_coupling_strength_median_rms_logregress(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_median_rms_logregress(:, ipoc))));

    sensitivity_coupling_strength_mean_MI(ipoc) = length(find(pvalues_coupling_strength_mean_MI(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_mean_MI(:, ipoc))));
    sensitivity_coupling_strength_mean_POS(ipoc) = length(find(pvalues_coupling_strength_mean_POS(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_mean_POS(:, ipoc))));
    sensitivity_coupling_strength_mean_z_rayleigh(ipoc) = length(find(pvalues_coupling_strength_mean_z_rayleigh(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_mean_z_rayleigh(:, ipoc))));
    sensitivity_coupling_strength_mean_U2watson(ipoc) = length(find(pvalues_coupling_strength_mean_U2watson(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_mean_U2watson(:, ipoc))));
    sensitivity_coupling_strength_mean_rms_logregress(ipoc) = length(find(pvalues_coupling_strength_mean_rms_logregress(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_coupling_strength_mean_rms_logregress(:, ipoc))));
    
end

%%% Show FP rate

figure; set(gcf, 'Color', 'w');
bar(1:5',[fp_coupling_strength_median_rms_logregress fp_coupling_strength_mean_rms_logregress; fp_coupling_strength_median_POS fp_coupling_strength_mean_POS; ...
    fp_coupling_strength_median_U2watson fp_coupling_strength_mean_U2watson; fp_coupling_strength_median_z_rayleigh fp_coupling_strength_mean_z_rayleigh;...
    fp_coupling_strength_median_MI fp_coupling_strength_mean_MI], 'grouped');
hline(0.05, 'k');
legend('median', 'mean');
xticklabels({'Log. regr.', 'POS', 'Watson', 'Rayleigh', 'MI'})

%%% Show sensitivity
% for one level of phase-outcome coupling
phase_outcome_coupling = 15;
ind_poc = find(cfg_simulations.strengths_phase_outcome_coupling==phase_outcome_coupling)-1;

figure; set(gcf, 'Color', 'w');
bar(1:5',[sensitivity_coupling_strength_median_rms_logregress(ind_poc) sensitivity_coupling_strength_mean_rms_logregress(ind_poc); sensitivity_coupling_strength_median_POS(ind_poc) sensitivity_coupling_strength_mean_POS(ind_poc); ...
    sensitivity_coupling_strength_median_U2watson(ind_poc) sensitivity_coupling_strength_mean_U2watson(ind_poc); sensitivity_coupling_strength_median_z_rayleigh(ind_poc) sensitivity_coupling_strength_mean_z_rayleigh(ind_poc);...
    sensitivity_coupling_strength_median_MI(ind_poc) sensitivity_coupling_strength_mean_MI(ind_poc)], 'grouped');
hline(0.05, 'k');
legend('median', 'mean');
xticklabels({'Log. regr.', 'POS', 'Watson', 'Rayleigh', 'MI'})

%%% Sensitivity vs. False Positive rate

figure('Color', 'w');
plot(fp_coupling_strength_median_rms_logregress*100, sensitivity_coupling_strength_median_rms_logregress(ind_poc)*100, 'or');
hold on;
plot(fp_coupling_strength_mean_rms_logregress*100, sensitivity_coupling_strength_mean_rms_logregress(ind_poc)*100, 'og');
hold on;
plot(fp_coupling_strength_median_POS*100, sensitivity_coupling_strength_median_POS(ind_poc)*100, 'sr');
hold on;
plot(fp_coupling_strength_mean_POS*100, sensitivity_coupling_strength_mean_POS(ind_poc)*100, 'sg');
hold on;
plot(fp_coupling_strength_median_U2watson*100, sensitivity_coupling_strength_median_U2watson(ind_poc)*100, 'dr');
hold on;
plot(fp_coupling_strength_mean_U2watson*100, sensitivity_coupling_strength_mean_U2watson(ind_poc)*100, 'dg');
hold on;
plot(fp_coupling_strength_median_z_rayleigh*100, sensitivity_coupling_strength_median_z_rayleigh(ind_poc)*100, '^r');
hold on;
plot(fp_coupling_strength_mean_z_rayleigh*100, sensitivity_coupling_strength_mean_z_rayleigh(ind_poc)*100, '^g');
hold on;
plot(fp_coupling_strength_median_MI*100, sensitivity_coupling_strength_median_MI(ind_poc)*100, 'vr');
hold on;
plot(fp_coupling_strength_mean_MI*100, sensitivity_coupling_strength_mean_MI(ind_poc)*100, 'vg');
ylim([0 100]);
xlabel('%False Positives', 'FontSize', 15);
ylabel('Sensitivity 20% phase-outcome coupling', 'FontSize', 15);

end
