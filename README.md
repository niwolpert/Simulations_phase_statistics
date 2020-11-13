Simulations phase statistics - Nicolai Wolpert (2020)
=======================

When using this function in any published study, please cite: Wolpert, N, & Tallon‐Baudry, C. Evaluation of different statistical procedures to estimate coupling between oscillatory phase and behavioral response (in preparation)

Copyright (C) 2020, Laboratoire de Neurosciences Cognitives, Nicolai Wolpert, Catherine Tallon-Baudry. Email: nicolaiwolpert@gmail.com

DISCLAIMER: This code is provided without explicit or implicit guarantee, and without any form of technical support. The code is not intended to be used for clinical purposes. The functions are free to use and can be redistributed, modified and adapted, under the terms of the CC BY-NC-SA version of creative commons license (see https://creativecommons.org/licenses/).

Purpose
-------------
The aim of these scripts is to compare different statistical methods to evaluate the coupling between oscillatory phase and a binary behavioral outcome. We create semi-artificial datasets mimicking real experiments, where we inject a statistical link between a simulated behavioral outcome and the phases of physiological oscillations from 30 participants. These two oscillations are the brain alpha rhythm (8-13 hz) and a very slow oscillation at ~0.05 Hz corresponding to the gastric rhythm (see 'alpha_phase_all_subjects' and 'EGG_phase_all_subjects' provided with the code). We systematically vary the strength of phase-outcome coupling, the coupling mode (1:1 to 4:1), the overall number of trials, the relative number of trials in the two outcome conditions, and evaluate different strategies to estimate phase-outcome coupling chance level, as well as significance at the individual or group level.
These scripts perform the computations of the simulations and save and vizualize the results.

Instructions
-----------------------

1. The functions require quite heavy computing and therefore use parallel computing to speed up calculations. They therefore require Matlab's Parallel Computing Toolbox to be installed: https://fr.mathworks.com/products/parallel-computing.html

2. This script makes use of some functions of the Circular Statistics Toolbox: https://github.com/anne-urai/Tools/tree/master/CircStat2012a

3. All the functions are called from 'SCRIPT_phase_stats_simulations_main'. Data sets with oscillatory phase time series are provided ('alpha_phase_all_subjects' and 'EGG_phase_all_subjects'). Download the data sets and specify the path to the the root directory and the Circular Statistics toolbox.

4. Specify the statistical & simulation parameters in the configuration structure 'cfg_simulations'

5. Load oscillatory phase time series 'alpha_phase_all_subjects' or 'EGG_phase_all_subjects'. These contain the time series of 30 participants (12-15 minute length)

6. Run simulations and evaulate results:
- Compare statistical tests (Phase Opposition Sum, Circular Logistic Regression, the Watson test, Modulation Index, and the Rayleigh test), and methods to assess significance at the group level (t-test on empirical vs. chance, surrogate average and three individual p-value combination methods)
  Here, different circular tests are compared in terms of sensitivty and False Positive rate as a function of the strength of the injected phase-outcome coupling. This section also compares different methods to estimate significance on the group level, and compares permutation statistics with tabulates statistics
- Number of trials in total:
  In this section, we compare the dependency of circular tests on the total number of trials
- Relative number of observations:
  Here, we systematically vary the relative number of observations in the two outcome conditions, and estimate sensitivity and False Positive rate of the five tests each each relative number of observations. We additionally investigate how a resampling procedure to control for an imbalanced amount of trials can recover sensitivity
- Optimal number of bins for Modulation Index (MI):
 Because MI depends on the logarithm of hit rate, it cannot be computed in the case that a given phase bin contains no hit. Here, estimate the relationship between the total number of trials (hits+misses), number of phase bins and the likelihood that MI cannot be computed.
  
