# Introduction
The code is for the paper 'Simulation-based, Finite-sample Inference for Privatized Data' by Jordan Awan and Zhanyu Wang. https://arxiv.org/abs/2303.05328v4

We denote Repro Sample method as Repro, and parametric bootstrap as PB.

# Run codes and draw figures
*For those experiments taking longer time to run, e.g., Figure 5, we split their procedure into two steps and provide intermediate results saved in the folder `linear_logistic/results/`. To save time, one can skip the first step and use the second step to draw the figures.

## Table 1: 95% confidence interval for Bernoulli
Please run `table1/Repro_Binomial_bounded.R` for the Repro result, and `table1/Awan_binomialDP.R` for the result by Awan and Slavkovic (2020).

## Figure 2: 95% confidence set for location-scale normal
1. Please run `figure2/Repro_Normal_DepthRegion_GDP.R` in R to obtain the four subfigures: `*_100_1_region.pdf` where `*` is `mahalanobis`, `halfspace`, `simplicial`, and `spatial`. 

2. We also provide intermediate boundary results saved in `csv` files, and you can run `figure2/Repro_Normal_DepthRegion_GDP_draw.R` to draw the figures.

## Figure 4: Poisson distribution
Please run `figure3/Repro_Poisson.R` in R to obtain the two subfigures: `CIPoisson.pdf` and `widthPoisson.pdf`. 

## Table 2: Location-scale normal (compared to PB)
Please run `table2/Repro_Normal_mu.R` and `table2/Repro_Normal_sigma.R` in R to obtain the results (saved in `csv` files) of Repro, and run `PB_Normal.R` for the results of PB.

## Table 3: Location-scale normal (compared to Karwa and Vadhan)
Please run `table3/Repro_Normal_mu_compareKarwaVadhan.R` for the results of repro, and `table3/Karwa_Vadhan.R` for the results by Karwa and Vadhan (2018).

## Figure 5: Compare Repro with PB for hypothesis testing on the slope in linear regression with different levels of the true slope
1. Please run `linear_logistic/Repro_LinearRegression_HypothesisTesting.R` and `linear_logistic/PB_LinearRegression_HypothesisTesting.R` to obtain the results which will be used in Figure 5 and 6 and Appendix Figure 1 and 2.

2. Please run the corresponding cells in `linear_logistic/Heatmaps.ipynb` and obtain `LR_HT_gdp=1.pdf`.

## Figure 6: Compare Repro with PB for hypothesis testing on the slope in linear regression with different clamping regions
1. (If you haven't done this for Figure 5:) Please run `linear_logistic/Repro_LinearRegression_HypothesisTesting.R` and `linear_logistic/PB_LinearRegression_HypothesisTesting.R` to obtain the results, which will be used in Figure 5 and 6 and Appendix Figure 1 and 2.

2. Please run the corresponding cells in `linear_logistic/Heatmaps.ipynb` and obtain `LR_HT_compareclamp.pdf`. 

## Figure 7: Compare Repro with DP-CI-ERM on the width and coverage for the confidence intervals of the coefficient in logistic regression
1. Please run `linear_logistic/Repro_Logistic.R` for the results of Repro, and `linear_logistic/ERM_Logistic.R` for the results of DP-CI-ERM. Note that the `linear_logistic/Repro_Logistic.R` may take longer time to run, and if it stops before it finishes, you can rerun it and it automatically continues to run since the intermediate results are saved in `linear_logistic/results/logistic/`. The final results will be saved in `linear_logistic/results/summary/`. 

2. Then run `linear_logistic/Logistic_comparison.R` to generate the left subfigure of Figure 7 which is named as `logistic_width_comparison.pdf`, and for the right subfigure, please run the corresponding cells in `linear_logistic/Heatmaps.ipynb` and obtain `logistic_coverage.pdf`.

## Appendix Table 1: 95% confidence intervals for clamped exponential distribution with Laplace noise
Please run `appendix_table1/Inversion_Exponential.R` in R to obtain the results for the inversion method, `appendix_table1/Repro_Exponential.R` for repro, and `appendix_table1/PB_Exponential.R` for PB.

## Appendix Figure 1 and 2: Use Repro and PB for hypothesis testing on linear regression with different privacy constraints
1. (If you haven't done this for Figure 5:) Please run `linear_logistic/Repro_LinearRegression_HypothesisTesting.R` and `linear_logistic/PB_LinearRegression_HypothesisTesting.R` to obtain the results, which will be used in Figure 5 and 6 and Appendix Figure 1 and 2.

2. Please run the corresponding cells in `linear_logistic/Heatmaps.ipynb` and obtain `Repro_LR_HT.pdf` and `PB_LR_HT.pdf`.

## Appendix Figure 3: Sensitivity space 
Please run `appendix_figure3/SensitivityHull.R` in R. It will generate `sensitivityHull.pdf`.

## Appendix Table 2: 95% confidence intervals for private Bernoullis with unknown n
Please run `appendix_table2/Repro_Binomial_Unbounded.R` in R for the results.

## Appendix Table 3: Privatized Mann–Whitney test (ep-DP)
Please run `appendix_table3/Repro.R` in R for the results of repro and `appendix_table3/Couch.R` for Couch et al. (2019). Please use `appendix_table3/Kazan.R` to tune the parameter for the best results by Kazan et al. (2023).

## Appendix Table 4: Privatized Mann–Whitney test (ep-GDP)
Please run `appendix_table4/Repro.R` in R for the results of repro and `appendix_table4/Couch.R` for Couch et al. (2019). Please use `appendix_table4/Kazan.R` to tune the parameter for the best results by Kazan et al. (2023).

## Appendix Table 5: Privatized Mann–Whitney test (change R)
Please run `appendix_table5/Repro_compare_R.R` in R for the results of repro.
