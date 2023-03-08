# Introduction

The code is for the paper 'Simulation-based Confidence Intervals and Hypothesis Tests for Privatized Data'. We denote Repro Sample method as Repro, and parametric bootstrap as PB.

# Run codes and draw figures
*For those experiments taking longer time to run, e.g., Figure 4, we split their procedure into two steps and provide intermediate results saved in the folder `./results/`. To save time, one can skip the first step and use the second step to draw the figures.

## Figure 2: 95% confidence set for location-scale normal
Please run `Repro_Normal_DepthRegion.R` in R to obtain the four subfigures: `repro_normal_*.pdf` where `*` is `mahalanobis`, `halfspace`, `simplicial`, and `spatial`. 

## Figure 3: Poisson distribution
Please run `Repro_Poisson.R` in R to obtain the two subfigures: `CIPoisson.pdf` and `widthPoisson.pdf`. 

## Table 1: Location-scale normal
Please run `Repro_Normal_mu.R` and `Repro_Normal_sigma.R` in R to obtain the results of Repro, and run `PB_Normal.R` for the results of PB.

## Figure 4: Compare Repro with PB for hypothesis testing on the slope in linear regression with different levels of the true slope
1. Please run `Repro_LinearRegression_HypothesisTesting.R` and `PB_LinearRegression_HypothesisTesting.R` to obtain the results which will be used in Figure 4 and 5 and Appendix Figure 1 and 2.

2. Please run the corresponding cells in `Heatmaps.ipynb` and obtain `LR_HT_gdp=1.pdf`.

## Figure 5: Compare Repro with PB for hypothesis testing on the slope in linear regression with different clamping regions
1. (If you haven't done this for Figure 4:) Please run `Repro_LinearRegression_HypothesisTesting.R` and `PB_LinearRegression_HypothesisTesting.R` to obtain the results, which will be used in Figure 4 and 5 and Appendix Figure 1 and 2.

2. Please run the corresponding cells in `Heatmaps.ipynb` and obtain `LR_HT_compareclamp.pdf`. 

## Figure 6: Compare Repro with DP-CI-ERM on the width and coverage for the confidence intervals of the coefficient in logistic regression
1. Please run `Repro_Logistic.R` for the results of Repro, and `ERM_Logistic.R` for the results of DP-CI-ERM. Note that the `Repro_Logistic.R` may take longer time to run, and if it stops before it finishes, you can rerun it and it automatically continues to run since the intermediate results are saved in `./results/logistic/`. The final results will be saved in `./results/summary/`. 

2. Then run `Logistic_comparison.R` to generate the left subfigure of Figure 6 which is named as `logistic_width_comparison.pdf`, and for the right subfigure, please run the corresponding cells in `Heatmaps.ipynb` and obtain `logistic_coverage.pdf`.

## Appendix Table 1: 95% confidence intervals for clamped exponential distribution with Laplace noise
Please run `Inversion_Exponential.R` in R to obtain the results for the inversion method, and run `PB_Exponential.R` for PB.

## Appendix Figure 1 and 2: Use Repro and PB for hypothesis testing on linear regression with different privacy constraints
1. (If you haven't done this for Figure 4:) Please run `Repro_LinearRegression_HypothesisTesting.R` and `PB_LinearRegression_HypothesisTesting.R` to obtain the results, which will be used in Figure 4 and 5 and Appendix Figure 1 and 2.

2. Please run the corresponding cells in `Heatmaps.ipynb` and obtain `Repro_LR_HT.pdf` and `PB_LR_HT.pdf`.

## Appendix Figure 3: Sensitivity space 
Please run `SensitivityHull.R` in R. It will generate `sensitivityHull.pdf`.

## Appendix Table 2: 95% confidence intervals for private Bernoullis with unknown n
Please run `Repro_Binomial_Unbounded.R` in R for the results.


