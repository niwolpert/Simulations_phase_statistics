function merge_results_per_phaseoutcome_coupling(cfg_simulations, root_dir, nexperiments, signific_thresh)
%MERGE_RESULTS_PER_PHASEOUTCOME_COUPLING
%   Merges results from the computations on varying strengths of
%   phase-outcome coupling.
%   - Computes p-values on the group level using different methods: 
%     Empirical vs. chance level, surrogate average, Fisher's p-value
%     combination, Stouffer's p-value combination and Edgington's p-value
%     combination
%   - Computes sensitivity for each statistical test and method to compute
%     group-level p-values as the proportion of experiments with a 
%     phase-outcome coupling higher than zero yielding a significant 
%     p-value (<.05).
%   - Computes False Positive rate for each statistical test and method to 
%     compute group-level p-values as the proportion of experiments 
%     yielding a significant p-value (<.05) in the absence of an injected
%     effect (phase-outcome coupling = 0).
%   
%   INPUTS
%   - cfg_simulations:      Configuration structure with simulation 
%                           parameters
%   - root_dir:                        Root directory where results will be
%                                      saved
%   - nexperiments:         Number of virtual experiments
%   - signific_thresh:      Threshold for significance (default: 0.05)
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

if nargin<4
    signific_thresh = 0.05;
end

fprintf('\n### Merging files across experiments, computing fp-rate and sensitivity...\n');

empirical_vs_chance_MI = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
empirical_vs_chance_POS = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
empirical_vs_chance_z_rayleigh = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
empirical_vs_chance_U2watson = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
empirical_vs_chance_rms_logregress = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));

pvalues_empirical_vs_chance_MI = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_empirical_vs_chance_POS = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_empirical_vs_chance_z_rayleigh = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_empirical_vs_chance_U2watson = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_empirical_vs_chance_rms_logregress = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));

pvalues_surravg_MI = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_surravg_POS = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_surravg_z_rayleigh = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_surravg_U2watson = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_surravg_rms_logregress = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));

pvalues_combined_stouffer_MI = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_stouffer_POS = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_stouffer_z_rayleigh = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_stouffer_U2watson = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_stouffer_rms_logregress = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));

pvalues_combined_fisher_MI = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_fisher_POS = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_fisher_z_rayleigh = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_fisher_U2watson = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_fisher_rms_logregress = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));

pvalues_combined_edington_MI = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_edington_POS = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_edington_z_rayleigh = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_edington_U2watson = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_edington_rms_logregress = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));

pvalues_combined_tippett_MI = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_tippett_POS = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_tippett_z_rayleigh = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_tippett_U2watson = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_tippett_rms_logregress = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));

pvalues_combined_friston_MI = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_friston_POS = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_friston_z_rayleigh = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_friston_U2watson = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));
pvalues_combined_friston_rms_logregress = nan(nexperiments, length(cfg_simulations.strengths_phase_outcome_coupling));

dispstat('','init');
for iexp=1:nexperiments
    
    dispstat(sprintf('Experiment %d%',iexp));
    
    if isfile(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_experiment' num2str(iexp) '.mat']))
        
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_empirical_experiment' num2str(iexp) '.mat']), 'MIs_empirical');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_empirical_experiment' num2str(iexp) '.mat']), 'POS_empirical');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_empirical_experiment' num2str(iexp) '.mat']), 'z_rayleigh_empirical');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_empirical_experiment' num2str(iexp) '.mat']), 'U2watson_empirical');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Empirical' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_empirical_experiment' num2str(iexp) '.mat']), 'rms_logregress_empirical');

        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_chance_levels_experiment' num2str(iexp) '.mat']), 'MIs_chance_levels');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_chance_levels_experiment' num2str(iexp) '.mat']), 'POS_chance_levels');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_chance_levels_experiment' num2str(iexp) '.mat']), 'z_rayleigh_chance_levels');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_chance_levels_experiment' num2str(iexp) '.mat']), 'U2watson_chance_levels');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Chance_levels' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_chance_levels_experiment' num2str(iexp) '.mat']), 'rms_logregress_chance_levels');

        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_distr_experiment' num2str(iexp) '.mat']), 'MIs_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_distr_experiment' num2str(iexp) '.mat']), 'POS_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_surr_distr_experiment' num2str(iexp) '.mat']), 'z_rayleigh_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_distr_experiment' num2str(iexp) '.mat']), 'U2watson_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_distr_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_distr');
        
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues_direct' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rayleigh_pvalues_direct_experiment' num2str(iexp) '.mat']), 'rayleigh_pvalues_direct');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues_direct' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_watson_pvalues_direct_experiment' num2str(iexp) '.mat']), 'watson_pvalues_direct');
        
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_MIs_surr_avg_distr_experiment' num2str(iexp) '.mat']), 'MIs_surr_avg_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_POS_surr_avg_distr_experiment' num2str(iexp) '.mat']), 'POS_surr_avg_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_z_rayleigh_surr_avg_distr_experiment' num2str(iexp) '.mat']), 'z_rayleigh_avg_surr_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_U2watson_surr_avg_distr_experiment' num2str(iexp) '.mat']), 'U2watson_surr_avg_distr');
        load(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Surrogate_avg_distributions' filesep cfg_simulations.oscillation filesep cfg_simulations.oscillation '_rms_logregress_surr_avg_distr_experiment' num2str(iexp) '.mat']), 'rms_logregress_surr_avg_distr');
        
        for ipoc=1:length(cfg_simulations.strengths_phase_outcome_coupling)
            
            empirical_vs_chance_MI(iexp, ipoc) = mean(MIs_empirical(:, ipoc))-mean(MIs_chance_levels(:, ipoc));
            empirical_vs_chance_POS(iexp, ipoc) = mean(POS_empirical(:, ipoc))-mean(POS_chance_levels(:, ipoc));
            empirical_vs_chance_z_rayleigh(iexp, ipoc) = mean(z_rayleigh_empirical(:, ipoc))-mean(z_rayleigh_chance_levels(:, ipoc));
            empirical_vs_chance_U2watson(iexp, ipoc) = mean(U2watson_empirical(:, ipoc))-mean(U2watson_chance_levels(:, ipoc));
            empirical_vs_chance_rms_logregress(iexp, ipoc) = mean(rms_logregress_empirical(:, ipoc))-mean(rms_logregress_chance_levels(:, ipoc));
            
            [~,pvalues_empirical_vs_chance_MI(iexp, ipoc)] = ttest(MIs_empirical(:, ipoc), MIs_chance_levels(:, ipoc),'Tail','right');
            [~,pvalues_empirical_vs_chance_POS(iexp, ipoc)] = ttest(POS_empirical(:, ipoc), POS_chance_levels(:, ipoc),'Tail','right');
            [~,pvalues_empirical_vs_chance_z_rayleigh(iexp, ipoc)] = ttest(z_rayleigh_empirical(:, ipoc), z_rayleigh_chance_levels(:, ipoc),'Tail','right');
            [~,pvalues_empirical_vs_chance_U2watson(iexp, ipoc)] = ttest(U2watson_empirical(:, ipoc), U2watson_chance_levels(:, ipoc),'Tail','right');
            [~,pvalues_empirical_vs_chance_rms_logregress(iexp, ipoc)] = ttest(rms_logregress_empirical(:, ipoc), rms_logregress_chance_levels(:, ipoc),'Tail','right');
            
            pvalues_surravg_MI(iexp, ipoc) = length(find(MIs_surr_avg_distr{ipoc}>mean(MIs_empirical(:, ipoc))))/length(MIs_surr_avg_distr{ipoc});
            pvalues_surravg_POS(iexp, ipoc) =  length(find(POS_surr_avg_distr{ipoc}>mean(POS_empirical(:, ipoc))))/length(POS_surr_avg_distr{ipoc});
            pvalues_surravg_z_rayleigh(iexp, ipoc) =  length(find(z_rayleigh_avg_surr_distr{ipoc}>mean(z_rayleigh_empirical(:, ipoc))))/length(z_rayleigh_avg_surr_distr{ipoc});
            pvalues_surravg_U2watson(iexp, ipoc) =  length(find(U2watson_surr_avg_distr{ipoc}>mean(U2watson_empirical(:, ipoc))))/length(U2watson_surr_avg_distr{ipoc});
            pvalues_surravg_rms_logregress(iexp, ipoc) =  length(find(rms_logregress_surr_avg_distr{ipoc}>mean(rms_logregress_empirical(:, ipoc))))/length(rms_logregress_surr_avg_distr{ipoc});
            
            nsubjects = size(MIs_empirical, 1);
            
            % compute combined p-values
            pvalues_individual_MI = nan(nsubjects, 1);
            pvalues_individual_POS = nan(nsubjects, 1);
            pvalues_individual_rayleigh = nan(nsubjects, 1);
            pvalues_individual_U2watson = nan(nsubjects, 1);
            pvalues_individual_rms_logregress = nan(nsubjects, 1);
            for isubject=1:nsubjects
                
                pvalues_individual_MI(isubject) = length(find(MIs_surr_distr{isubject, ipoc}>MIs_empirical(isubject, ipoc)))/length(MIs_surr_distr{isubject, ipoc});
                if pvalues_individual_MI(isubject)<(1/(2*cfg_simulations.nperm))
                    pvalues_individual_MI(isubject) = (1/(2*cfg_simulations.nperm));
                end
                pvalues_individual_POS(isubject) = length(find(POS_surr_distr{isubject, ipoc}>POS_empirical(isubject, ipoc)))/length(POS_surr_distr{isubject, ipoc});
                if pvalues_individual_POS(isubject)<(1/(2*cfg_simulations.nperm))
                    pvalues_individual_POS(isubject) = (1/(2*cfg_simulations.nperm));
                end
                pvalues_individual_rayleigh(isubject) = length(find(z_rayleigh_surr_distr{isubject, ipoc}>z_rayleigh_empirical(isubject, ipoc)))/length(z_rayleigh_surr_distr{isubject, ipoc});
                if pvalues_individual_rayleigh(isubject)<(1/(2*cfg_simulations.nperm))
                    pvalues_individual_rayleigh(isubject) = (1/(2*cfg_simulations.nperm));
                end
                pvalues_individual_U2watson(isubject) = length(find(U2watson_surr_distr{isubject, ipoc}>U2watson_empirical(isubject, ipoc)))/length(U2watson_surr_distr{isubject, ipoc});
                if pvalues_individual_U2watson(isubject)<(1/(2*cfg_simulations.nperm))
                    pvalues_individual_U2watson(isubject) = (1/(2*cfg_simulations.nperm));
                end
                pvalues_individual_rms_logregress(isubject) = length(find(rms_logregress_surr_distr{isubject, ipoc}>rms_logregress_empirical(isubject, ipoc)))/length(rms_logregress_surr_distr{isubject, ipoc});
                if pvalues_individual_rms_logregress(isubject)<(1/(2*cfg_simulations.nperm))
                    pvalues_individual_rms_logregress(isubject) = (1/(2*cfg_simulations.nperm));
                end
            end
            
            pvalues_combined_stouffer_MI(iexp, ipoc) = squeeze(1-normcdf(sum(norminv(1-pvalues_individual_MI),1)./sqrt(size(pvalues_individual_MI,1))));
            pvalues_combined_stouffer_POS(iexp, ipoc) = squeeze(1-normcdf(sum(norminv(1-pvalues_individual_POS),1)./sqrt(size(pvalues_individual_POS,1))));
            pvalues_combined_stouffer_z_rayleigh(iexp, ipoc) = squeeze(1-normcdf(sum(norminv(1-pvalues_individual_rayleigh),1)./sqrt(size(pvalues_individual_rayleigh,1))));
            pvalues_combined_stouffer_U2watson(iexp, ipoc) = squeeze(1-normcdf(sum(norminv(1-pvalues_individual_U2watson),1)./sqrt(size(pvalues_individual_U2watson,1))));
            pvalues_combined_stouffer_rms_logregress(iexp, ipoc) = squeeze(1-normcdf(sum(norminv(1-pvalues_individual_rms_logregress),1)./sqrt(size(pvalues_individual_rms_logregress,1))));

            pvalues_combined_fisher_MI(iexp, ipoc) = chi2cdf(squeeze(-2*sum(log(pvalues_individual_MI),1)),2*size(pvalues_individual_MI,1),'upper');
            pvalues_combined_fisher_POS(iexp, ipoc) = chi2cdf(squeeze(-2*sum(log(pvalues_individual_POS),1)),2*size(pvalues_individual_POS,1),'upper');
            pvalues_combined_fisher_z_rayleigh(iexp, ipoc) = chi2cdf(squeeze(-2*sum(log(pvalues_individual_rayleigh),1)),2*size(pvalues_individual_rayleigh,1),'upper');
            pvalues_combined_fisher_U2watson(iexp, ipoc) = chi2cdf(squeeze(-2*sum(log(pvalues_individual_U2watson),1)),2*size(pvalues_individual_U2watson,1),'upper');
            pvalues_combined_fisher_rms_logregress(iexp, ipoc) = chi2cdf(squeeze(-2*sum(log(pvalues_individual_rms_logregress),1)),2*size(pvalues_individual_rms_logregress,1),'upper');

            pvalues_combined_edington_MI(iexp, ipoc) = squeeze(normcdf(sum(pvalues_individual_MI,1),0.5*size(pvalues_individual_MI,1),sqrt(size(pvalues_individual_MI,1)/12)));
            pvalues_combined_edington_POS(iexp, ipoc) = squeeze(normcdf(sum(pvalues_individual_POS,1),0.5*size(pvalues_individual_POS,1),sqrt(size(pvalues_individual_POS,1)/12)));
            pvalues_combined_edington_z_rayleigh(iexp, ipoc) = squeeze(normcdf(sum(pvalues_individual_rayleigh,1),0.5*size(pvalues_individual_rayleigh,1),sqrt(size(pvalues_individual_rayleigh,1)/12)));
            pvalues_combined_edington_U2watson(iexp, ipoc) = squeeze(normcdf(sum(pvalues_individual_U2watson,1),0.5*size(pvalues_individual_U2watson,1),sqrt(size(pvalues_individual_U2watson,1)/12)));
            pvalues_combined_edington_rms_logregress(iexp, ipoc) = squeeze(normcdf(sum(pvalues_individual_rms_logregress,1),0.5*size(pvalues_individual_rms_logregress,1),sqrt(size(pvalues_individual_rms_logregress,1)/12)));

            pvalues_combined_tippett_MI(iexp, ipoc) = squeeze(1-(1-min(pvalues_individual_MI,[],1)).^size(pvalues_individual_MI,1));
            pvalues_combined_tippett_POS(iexp, ipoc) = squeeze(1-(1-min(pvalues_individual_POS,[],1)).^size(pvalues_individual_POS,1));
            pvalues_combined_tippett_z_rayleigh(iexp, ipoc) = squeeze(1-(1-min(pvalues_individual_rayleigh,[],1)).^size(pvalues_individual_rayleigh,1));
            pvalues_combined_tippett_U2watson(iexp, ipoc) = squeeze(1-(1-min(pvalues_individual_U2watson,[],1)).^size(pvalues_individual_U2watson,1));
            pvalues_combined_tippett_rms_logregress(iexp, ipoc) = squeeze(1-(1-min(pvalues_individual_rms_logregress,[],1)).^size(pvalues_individual_rms_logregress,1));

            pvalues_combined_friston_MI(iexp, ipoc) = squeeze((max(pvalues_individual_MI,[],1)).^size(pvalues_individual_MI,1));
            pvalues_combined_friston_POS(iexp, ipoc) = squeeze((max(pvalues_individual_POS,[],1)).^size(pvalues_individual_POS,1));
            pvalues_combined_friston_z_rayleigh(iexp, ipoc) = squeeze((max(pvalues_individual_rayleigh,[],1)).^size(pvalues_individual_rayleigh,1));
            pvalues_combined_friston_U2watson(iexp, ipoc) = squeeze((max(pvalues_individual_U2watson,[],1)).^size(pvalues_individual_U2watson,1));
            pvalues_combined_friston_rms_logregress(iexp, ipoc) = squeeze((max(pvalues_individual_rms_logregress,[],1)).^size(pvalues_individual_rms_logregress,1));
            
        end
    end
end

% compute false positive (FP) rate

fp_empirical_vs_chance_MI = length(find(pvalues_empirical_vs_chance_MI(:, 1)<signific_thresh))/length(find(~isnan(pvalues_empirical_vs_chance_MI(:, 1))));
fp_empirical_vs_chance_POS = length(find(pvalues_empirical_vs_chance_POS(:, 1)<signific_thresh))/length(find(~isnan(pvalues_empirical_vs_chance_POS(:, 1))));
fp_empirical_vs_chance_z_rayleigh = length(find(pvalues_empirical_vs_chance_z_rayleigh(:, 1)<signific_thresh))/length(find(~isnan(pvalues_empirical_vs_chance_z_rayleigh(:, 1))));
fp_empirical_vs_chance_U2watson = length(find(pvalues_empirical_vs_chance_U2watson(:, 1)<signific_thresh))/length(find(~isnan(pvalues_empirical_vs_chance_U2watson(:, 1))));
fp_empirical_vs_chance_rms_logregress = length(find(pvalues_empirical_vs_chance_rms_logregress(:, 1)<signific_thresh))/length(find(~isnan(pvalues_empirical_vs_chance_rms_logregress(:, 1))));

fp_surravg_MI = length(find(pvalues_surravg_MI(:, 1)<signific_thresh))/length(find(~isnan(pvalues_surravg_MI(:, 1))));
fp_surravg_POS = length(find(pvalues_surravg_POS(:, 1)<signific_thresh))/length(find(~isnan(pvalues_surravg_POS(:, 1))));
fp_surravg_z_rayleigh = length(find(pvalues_surravg_z_rayleigh(:, 1)<signific_thresh))/length(find(~isnan(pvalues_surravg_z_rayleigh(:, 1))));
fp_surravg_U2watson = length(find(pvalues_surravg_U2watson(:, 1)<signific_thresh))/length(find(~isnan(pvalues_surravg_U2watson(:, 1))));
fp_surravg_rms_logregress = length(find(pvalues_surravg_rms_logregress(:, 1)<signific_thresh))/length(find(~isnan(pvalues_surravg_rms_logregress(:, 1))));

fp_combined_stouffer_MI = length(find(pvalues_combined_stouffer_MI(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_stouffer_MI(:, 1))));
fp_combined_stouffer_POS = length(find(pvalues_combined_stouffer_POS(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_stouffer_POS(:, 1))));
fp_combined_stouffer_z_rayleigh = length(find(pvalues_combined_stouffer_z_rayleigh(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_stouffer_z_rayleigh(:, 1))));
fp_combined_stouffer_U2watson = length(find(pvalues_combined_stouffer_U2watson(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_stouffer_U2watson(:, 1))));
fp_combined_stouffer_rms_logregress = length(find(pvalues_combined_stouffer_rms_logregress(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_stouffer_rms_logregress(:, 1))));

fp_combined_fisher_MI = length(find(pvalues_combined_fisher_MI(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_fisher_MI(:, 1))));
fp_combined_fisher_POS = length(find(pvalues_combined_fisher_POS(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_fisher_POS(:, 1))));
fp_combined_fisher_z_rayleigh = length(find(pvalues_combined_fisher_z_rayleigh(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_fisher_z_rayleigh(:, 1))));
fp_combined_fisher_U2watson = length(find(pvalues_combined_fisher_U2watson(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_fisher_U2watson(:, 1))));
fp_combined_fisher_rms_logregress = length(find(pvalues_combined_fisher_rms_logregress(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_fisher_rms_logregress(:, 1))));

fp_combined_edington_MI = length(find(pvalues_combined_edington_MI(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_edington_MI(:, 1))));
fp_combined_edington_POS = length(find(pvalues_combined_edington_POS(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_edington_POS(:, 1))));
fp_combined_edington_z_rayleigh = length(find(pvalues_combined_edington_z_rayleigh(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_edington_z_rayleigh(:, 1))));
fp_combined_edington_U2watson = length(find(pvalues_combined_edington_U2watson(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_edington_U2watson(:, 1))));
fp_combined_edington_rms_logregress = length(find(pvalues_combined_edington_rms_logregress(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_edington_rms_logregress(:, 1))));

fp_combined_tippett_MI = length(find(pvalues_combined_tippett_MI(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_tippett_MI(:, 1))));
fp_combined_tippett_POS = length(find(pvalues_combined_tippett_POS(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_tippett_POS(:, 1))));
fp_combined_tippett_z_rayleigh = length(find(pvalues_combined_tippett_z_rayleigh(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_tippett_z_rayleigh(:, 1))));
fp_combined_tippett_U2watson = length(find(pvalues_combined_tippett_U2watson(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_tippett_U2watson(:, 1))));
fp_combined_tippett_rms_logregress = length(find(pvalues_combined_tippett_rms_logregress(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_tippett_rms_logregress(:, 1))));

fp_combined_friston_MI = length(find(pvalues_combined_friston_MI(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_friston_MI(:, 1))));
fp_combined_friston_POS = length(find(pvalues_combined_friston_POS(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_friston_POS(:, 1))));
fp_combined_friston_z_rayleigh = length(find(pvalues_combined_friston_z_rayleigh(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_friston_z_rayleigh(:, 1))));
fp_combined_friston_U2watson = length(find(pvalues_combined_friston_U2watson(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_friston_U2watson(:, 1))));
fp_combined_friston_rms_logregress = length(find(pvalues_combined_friston_rms_logregress(:, 1)<signific_thresh))/length(find(~isnan(pvalues_combined_friston_rms_logregress(:, 1))));

% compute sensitivity
sensitivity_empirical_vs_chance_MI = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_empirical_vs_chance_POS = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_empirical_vs_chance_z_rayleigh = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_empirical_vs_chance_U2watson = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_empirical_vs_chance_rms_logregress = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));

sensitivity_surravg_MI = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_surravg_POS = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_surravg_z_rayleigh = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_surravg_U2watson = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_surravg_rms_logregress = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));

sensitivity_combined_stouffer_MI = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_stouffer_POS = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_stouffer_z_rayleigh = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_stouffer_U2watson = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_stouffer_rms_logregress = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));

sensitivity_combined_fisher_MI = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_fisher_POS = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_fisher_z_rayleigh = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_fisher_U2watson = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_fisher_rms_logregress = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));

sensitivity_combined_edington_MI = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_edington_POS = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_edington_z_rayleigh = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_edington_U2watson = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_edington_rms_logregress = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));

sensitivity_combined_tippett_MI = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_tippett_POS = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_tippett_z_rayleigh = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_tippett_U2watson = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_tippett_rms_logregress = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));

sensitivity_combined_friston_MI = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_friston_POS = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_friston_z_rayleigh = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_friston_U2watson = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));
sensitivity_combined_friston_rms_logregress = nan(1, length(cfg_simulations.strengths_phase_outcome_coupling));

for ipoc=1:length(cfg_simulations.strengths_phase_outcome_coupling)
    
    sensitivity_empirical_vs_chance_MI(ipoc) = length(find(pvalues_empirical_vs_chance_MI(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_empirical_vs_chance_MI(:, ipoc))));
    sensitivity_empirical_vs_chance_POS(ipoc) = length(find(pvalues_empirical_vs_chance_POS(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_empirical_vs_chance_POS(:, ipoc))));
    sensitivity_empirical_vs_chance_z_rayleigh(ipoc) = length(find(pvalues_empirical_vs_chance_z_rayleigh(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_empirical_vs_chance_z_rayleigh(:, ipoc))));
    sensitivity_empirical_vs_chance_U2watson(ipoc) = length(find(pvalues_empirical_vs_chance_U2watson(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_empirical_vs_chance_U2watson(:, ipoc))));
    sensitivity_empirical_vs_chance_rms_logregress(ipoc) = length(find(pvalues_empirical_vs_chance_rms_logregress(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_empirical_vs_chance_rms_logregress(:, ipoc))));
    
    sensitivity_surravg_MI(ipoc) = length(find(pvalues_surravg_MI(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_surravg_MI(:, ipoc))));
    sensitivity_surravg_POS(ipoc) = length(find(pvalues_surravg_POS(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_surravg_POS(:, ipoc))));
    sensitivity_surravg_z_rayleigh(ipoc) = length(find(pvalues_surravg_z_rayleigh(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_surravg_z_rayleigh(:, ipoc))));
    sensitivity_surravg_U2watson(ipoc) = length(find(pvalues_surravg_U2watson(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_surravg_U2watson(:, ipoc))));
    sensitivity_surravg_rms_logregress(ipoc) = length(find(pvalues_surravg_rms_logregress(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_surravg_rms_logregress(:, ipoc))));
    
    sensitivity_combined_stouffer_MI(ipoc) = length(find(pvalues_combined_stouffer_MI(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_stouffer_MI(:, ipoc))));
    sensitivity_combined_stouffer_POS(ipoc) = length(find(pvalues_combined_stouffer_POS(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_stouffer_POS(:, ipoc))));
    sensitivity_combined_stouffer_z_rayleigh(ipoc) = length(find(pvalues_combined_stouffer_z_rayleigh(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_stouffer_z_rayleigh(:, ipoc))));
    sensitivity_combined_stouffer_U2watson(ipoc) = length(find(pvalues_combined_stouffer_U2watson(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_stouffer_U2watson(:, ipoc))));
    sensitivity_combined_stouffer_rms_logregress(ipoc) = length(find(pvalues_combined_stouffer_rms_logregress(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_stouffer_rms_logregress(:, ipoc))));
    
    sensitivity_combined_fisher_MI(ipoc) = length(find(pvalues_combined_fisher_MI(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_fisher_MI(:, ipoc))));
    sensitivity_combined_fisher_POS(ipoc) = length(find(pvalues_combined_fisher_POS(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_fisher_POS(:, ipoc))));
    sensitivity_combined_fisher_z_rayleigh(ipoc) = length(find(pvalues_combined_fisher_z_rayleigh(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_fisher_z_rayleigh(:, ipoc))));
    sensitivity_combined_fisher_U2watson(ipoc) = length(find(pvalues_combined_fisher_U2watson(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_fisher_U2watson(:, ipoc))));
    sensitivity_combined_fisher_rms_logregress(ipoc) = length(find(pvalues_combined_fisher_rms_logregress(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_fisher_rms_logregress(:, ipoc))));
    
    sensitivity_combined_edington_MI(ipoc) = length(find(pvalues_combined_edington_MI(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_edington_MI(:, ipoc))));
    sensitivity_combined_edington_POS(ipoc) = length(find(pvalues_combined_edington_POS(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_edington_POS(:, ipoc))));
    sensitivity_combined_edington_z_rayleigh(ipoc) = length(find(pvalues_combined_edington_z_rayleigh(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_edington_z_rayleigh(:, ipoc))));
    sensitivity_combined_edington_U2watson(ipoc) = length(find(pvalues_combined_edington_U2watson(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_edington_U2watson(:, ipoc))));
    sensitivity_combined_edington_rms_logregress(ipoc) = length(find(pvalues_combined_edington_rms_logregress(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_edington_rms_logregress(:, ipoc))));
    
    sensitivity_combined_tippett_MI(ipoc) = length(find(pvalues_combined_tippett_MI(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_tippett_MI(:, ipoc))));
    sensitivity_combined_tippett_POS(ipoc) = length(find(pvalues_combined_tippett_POS(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_tippett_POS(:, ipoc))));
    sensitivity_combined_tippett_z_rayleigh(ipoc) = length(find(pvalues_combined_tippett_z_rayleigh(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_tippett_z_rayleigh(:, ipoc))));
    sensitivity_combined_tippett_U2watson(ipoc) = length(find(pvalues_combined_tippett_U2watson(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_tippett_U2watson(:, ipoc))));
    sensitivity_combined_tippett_rms_logregress(ipoc) = length(find(pvalues_combined_tippett_rms_logregress(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_tippett_rms_logregress(:, ipoc))));
    
    sensitivity_combined_friston_MI(ipoc) = length(find(pvalues_combined_friston_MI(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_friston_MI(:, ipoc))));
    sensitivity_combined_friston_POS(ipoc) = length(find(pvalues_combined_friston_POS(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_friston_POS(:, ipoc))));
    sensitivity_combined_friston_z_rayleigh(ipoc) = length(find(pvalues_combined_friston_z_rayleigh(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_friston_z_rayleigh(:, ipoc))));
    sensitivity_combined_friston_U2watson(ipoc) = length(find(pvalues_combined_friston_U2watson(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_friston_U2watson(:, ipoc))));
    sensitivity_combined_friston_rms_logregress(ipoc) = length(find(pvalues_combined_friston_rms_logregress(:, ipoc)<signific_thresh))/length(find(~isnan(pvalues_combined_friston_rms_logregress(:, ipoc))));
    
end

%%% save

% p-values

if ~isdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation filesep]))
end

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_empirical_vs_chance_MI.mat']), 'pvalues_empirical_vs_chance_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_empirical_vs_chance_POS.mat']), 'pvalues_empirical_vs_chance_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_empirical_vs_chance_z_rayleigh.mat']), 'pvalues_empirical_vs_chance_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_empirical_vs_chance_U2watson.mat']), 'pvalues_empirical_vs_chance_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_empirical_vs_chance_rms_logregress.mat']), 'pvalues_empirical_vs_chance_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_surravg_MI.mat']), 'pvalues_surravg_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_surravg_POS.mat']), 'pvalues_surravg_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_surravg_z_rayleigh.mat']), 'pvalues_surravg_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_surravg_U2watson.mat']), 'pvalues_surravg_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_surravg_rms_logregress.mat']), 'pvalues_surravg_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_stouffer_MI.mat']), 'pvalues_combined_stouffer_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_stouffer_POS.mat']), 'pvalues_combined_stouffer_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_stouffer_z_rayleigh.mat']), 'pvalues_combined_stouffer_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_stouffer_U2watson.mat']), 'pvalues_combined_stouffer_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_stouffer_rms_logregress.mat']), 'pvalues_combined_stouffer_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_fisher_MI.mat']), 'pvalues_combined_fisher_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_fisher_POS.mat']), 'pvalues_combined_fisher_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_fisher_z_rayleigh.mat']), 'pvalues_combined_fisher_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_fisher_U2watson.mat']), 'pvalues_combined_fisher_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_fisher_rms_logregress.mat']), 'pvalues_combined_fisher_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_edington_MI.mat']), 'pvalues_combined_edington_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_edington_POS.mat']), 'pvalues_combined_edington_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_edington_z_rayleigh.mat']), 'pvalues_combined_edington_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_edington_U2watson.mat']), 'pvalues_combined_edington_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_edington_rms_logregress.mat']), 'pvalues_combined_edington_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_tippett_MI.mat']), 'pvalues_combined_tippett_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_tippett_POS.mat']), 'pvalues_combined_tippett_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_tippett_z_rayleigh.mat']), 'pvalues_combined_tippett_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_tippett_U2watson.mat']), 'pvalues_combined_tippett_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_tippett_rms_logregress.mat']), 'pvalues_combined_tippett_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_friston_MI.mat']), 'pvalues_combined_friston_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_friston_POS.mat']), 'pvalues_combined_friston_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_friston_z_rayleigh.mat']), 'pvalues_combined_friston_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_friston_U2watson.mat']), 'pvalues_combined_friston_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Pvalues' filesep cfg_simulations.oscillation '_pvalues_combined_friston_rms_logregress.mat']), 'pvalues_combined_friston_rms_logregress');


% FP-rate

if ~isdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation filesep]))
end

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_MI.mat']), 'fp_empirical_vs_chance_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_POS.mat']), 'fp_empirical_vs_chance_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_z_rayleigh.mat']), 'fp_empirical_vs_chance_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_U2watson.mat']), 'fp_empirical_vs_chance_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_empirical_vs_chance_rms_logregress.mat']), 'fp_empirical_vs_chance_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_MI.mat']), 'fp_surravg_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_POS.mat']), 'fp_surravg_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_z_rayleigh.mat']), 'fp_surravg_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_U2watson.mat']), 'fp_surravg_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_surravg_rms_logregress.mat']), 'fp_surravg_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_MI.mat']), 'fp_combined_stouffer_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_POS.mat']), 'fp_combined_stouffer_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_z_rayleigh.mat']), 'fp_combined_stouffer_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_U2watson.mat']), 'fp_combined_stouffer_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_stouffer_rms_logregress.mat']), 'fp_combined_stouffer_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_MI.mat']), 'fp_combined_fisher_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_POS.mat']), 'fp_combined_fisher_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_z_rayleigh.mat']), 'fp_combined_fisher_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_U2watson.mat']), 'fp_combined_fisher_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_fisher_rms_logregress.mat']), 'fp_combined_fisher_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_MI.mat']), 'fp_combined_edington_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_POS.mat']), 'fp_combined_edington_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_z_rayleigh.mat']), 'fp_combined_edington_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_U2watson.mat']), 'fp_combined_edington_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_edington_rms_logregress.mat']), 'fp_combined_edington_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_MI.mat']), 'fp_combined_tippett_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_POS.mat']), 'fp_combined_tippett_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_z_rayleigh.mat']), 'fp_combined_tippett_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_U2watson.mat']), 'fp_combined_tippett_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_tippett_rms_logregress.mat']), 'fp_combined_tippett_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_MI.mat']), 'fp_combined_friston_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_POS.mat']), 'fp_combined_friston_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_z_rayleigh.mat']), 'fp_combined_friston_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_U2watson.mat']), 'fp_combined_friston_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'FPrate' filesep cfg_simulations.oscillation '_fp_combined_friston_rms_logregress.mat']), 'fp_combined_friston_rms_logregress');

% Sensitivity

if ~isdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep]))
    mkdir(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation filesep]))
end

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_MI.mat']), 'sensitivity_empirical_vs_chance_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_POS.mat']), 'sensitivity_empirical_vs_chance_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_z_rayleigh.mat']), 'sensitivity_empirical_vs_chance_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_U2watson.mat']), 'sensitivity_empirical_vs_chance_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_empirical_vs_chance_rms_logregress.mat']), 'sensitivity_empirical_vs_chance_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_MI.mat']), 'sensitivity_surravg_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_POS.mat']), 'sensitivity_surravg_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_z_rayleigh.mat']), 'sensitivity_surravg_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_U2watson.mat']), 'sensitivity_surravg_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_surravg_rms_logregress.mat']), 'sensitivity_surravg_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_MI.mat']), 'sensitivity_combined_stouffer_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_POS.mat']), 'sensitivity_combined_stouffer_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_z_rayleigh.mat']), 'sensitivity_combined_stouffer_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_U2watson.mat']), 'sensitivity_combined_stouffer_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_stouffer_rms_logregress.mat']), 'sensitivity_combined_stouffer_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_MI.mat']), 'sensitivity_combined_fisher_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_POS.mat']), 'sensitivity_combined_fisher_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_z_rayleigh.mat']), 'sensitivity_combined_fisher_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_U2watson.mat']), 'sensitivity_combined_fisher_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_fisher_rms_logregress.mat']), 'sensitivity_combined_fisher_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_MI.mat']), 'sensitivity_combined_edington_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_POS.mat']), 'sensitivity_combined_edington_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_z_rayleigh.mat']), 'sensitivity_combined_edington_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_U2watson.mat']), 'sensitivity_combined_edington_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_edington_rms_logregress.mat']), 'sensitivity_combined_edington_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_MI.mat']), 'sensitivity_combined_tippett_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_POS.mat']), 'sensitivity_combined_tippett_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_z_rayleigh.mat']), 'sensitivity_combined_tippett_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_U2watson.mat']), 'sensitivity_combined_tippett_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_tippett_rms_logregress.mat']), 'sensitivity_combined_tippett_rms_logregress');

save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_MI.mat']), 'sensitivity_combined_friston_MI');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_POS.mat']), 'sensitivity_combined_friston_POS');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_z_rayleigh.mat']), 'sensitivity_combined_friston_z_rayleigh');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_U2watson.mat']), 'sensitivity_combined_friston_U2watson');
save(strcat([root_dir 'Compare_sensitivity_methods' filesep 'Mode' num2str(cfg_simulations.coupling_mode) filesep 'Sensitivity' filesep cfg_simulations.oscillation '_sensitivity_combined_friston_rms_logregress.mat']), 'sensitivity_combined_friston_rms_logregress');

end