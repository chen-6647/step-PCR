##################

A stepwise principal component regression method for localizing effects and individual differences in hierarchical structures

##################


We provide R code for Simulations 1-3 and the empirical application.

##### Simulation 1
1. The file "data_and_models.R" generates simulated data and performs multiple regression modeling.
2. The file "analysis_p.R" aggregates t-tests and p-values from 1000 simulations.

##### Simulation 2
1. Each folder is named as "Condition [condition number]_[sample size]".
2. The file "data_and_models.R" generates simulated data and performs step-PCR modeling.
3. The file "analysis_model_comparison.R" compares the performances of different model-selection methods.

##### Simulation 3
1. The file "sim_parameters.R" generates a large pool of parameters used for simulations.
2. The files "Condition [condition number].R" generate simulated data and perform step-PCR modeling.
3. The files "results_[condition number].R" analyze and plot the results.

##### Empirical Application
1. The file "analysis.R" analyzes dark triad and hypersensitive narcissism data based on step-PCR.
2. Data is available at https://openpsychometrics.org/_rawdata/.

