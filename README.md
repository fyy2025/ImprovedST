# Computing Codes for the Paper "Improving the within-Node Estimation of Survival Trees while Retaining Interpretability"
### Haolin Li*, Yiyang Fan*, Jianwen Cai#


## Description

This repository contains computing codes for the paper "Improving the within-Node Estimation of Survival Trees while Retaining Interpretability". Please click [here](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13821) for the full text of the paper.

## Folders

### 1-Data Generation

In this folder, we summarize the computing codes for generating survival data following a Weibull hazard model. There are two subfolders, `Continuous` containing codes for generating data with only continuous variables, and `Categorcial` containing codes for data involving categorical variables. Each of the two folders contains two folders, `low dimension` where each dataset generated has 10 covariates and `high dimension` where each dataset generated has (n+50) covariates (n is the number of observations). In each folder, the names of the code and the corresponding simulation scenarios in the paper are as follows,


* *data_generation_1.r* - linear effect, 100 sample size, 80% event rate.
* *data_generation_2.r* - linear effect, 100 sample size, 50% event rate.
* *data_generation_3.r* - linear effect, 300 sample size, 50% event rate.
* *data_generation_4.r* - linear effect, 300 sample size, 30% event rate.
* *data_generation_5.r* - linear effect, 500 sample size, 30% event rate.
* *data_generation_6.r* - linear effect, 500 sample size, 20% event rate.
* *data_generation_7.r* - nonlinear covariate effect, 100 sample size, 80% event rate.
* *data_generation_8.r* - nonlinear covariate effect, 100 sample size, 50% event rate.
* *data_generation_9.r* - nonlinear covariate effect, 300 sample size, 50% event rate.
* *data_generation_10.r* - nonlinear covariate effect, 300 sample size, 30% event rate.
* *data_generation_11.r* - nonlinear covariate effect, 500 sample size, 30% event rate.
* *data_generation_12.r* - nonlinear covariate effect, 500 sample size, 20% event rate.
* *data_generation_13.r* - tree-based covariate effect, 100 sample size, 80% event rate.
* *data_generation_14.r* - tree-based covariate effect, 100 sample size, 50% event rate.
* *data_generation_15.r* - tree-based covariate effect, 300 sample size, 50% event rate.
* *data_generation_16.r* - tree-based covariate effect, 300 sample size, 30% event rate.
* *data_generation_17.r* - tree-based covariate effect, 500 sample size, 30% event rate.
* *data_generation_18.r* - tree-based covariate effect, 500 sample size, 20% event rate.

### 2-Analysis

In this folder, we summarize the computing codes for fitting the proposed survival tree ensemble method and comparing the concordance index with single tree approach. The structure is the same as that in `Data Generation` In each folder, the names of the code and the corresponding simulation scenarios in the paper are as follows,

* *analysis_1.r* - linear effect, 100 sample size, 80% event rate.
* *analysis_2.r* - linear effect, 100 sample size, 50% event rate.
* *analysis_3.r* - linear effect, 300 sample size, 50% event rate.
* *analysis_4.r* - linear effect, 300 sample size, 30% event rate.
* *analysis_5.r* - linear effect, 500 sample size, 30% event rate.
* *analysis_6.r* - linear effect, 500 sample size, 20% event rate.
* *analysis_7.r* - nonlinear covariate effect, 100 sample size, 80% event rate.
* *analysis_8.r* - nonlinear covariate effect, 100 sample size, 50% event rate.
* *analysis_9.r* - nonlinear covariate effect, 300 sample size, 50% event rate.
* *analysis_10.r* - nonlinear covariate effect, 300 sample size, 30% event rate.
* *analysis_11.r* - nonlinear covariate effect, 500 sample size, 30% event rate.
* *analysis_12.r* - nonlinear covariate effect, 500 sample size, 20% event rate.
* *analysis_13.r* - tree-based covariate effect, 100 sample size, 80% event rate.
* *analysis_14.r* - tree-based covariate effect, 100 sample size, 50% event rate.
* *analysis_15.r* - tree-based covariate effect, 300 sample size, 50% event rate.
* *analysis_16.r* - tree-based covariate effect, 300 sample size, 30% event rate.
* *analysis_17.r* - tree-based covariate effect, 500 sample size, 30% event rate.
* *analysis_18.r* - tree-based covariate effect, 500 sample size, 20% event rate.

### Data availability statement

All data are publicly available. The dataset used in Section 4.1 in the paper is included in the ‘survival’ package in R. The dataset from Section 4.2 is available under the Creative Commons Attribution 4.0 International (CC BY 4.0) license at https://plos.figshare.com/articles/Survival_analysis_of_heart_failure_patients_A_case_study/5227684/1. The dataset in Section 4.3 is available in the ‘curatedOvarianData’ package in R.

## References

* Li, H., Fan, Y., & Cai, J. (2024+). Improving the within-Node Estimation of Survival Trees while Retaining Interpretability. Manuscript Submitted for Publication.
