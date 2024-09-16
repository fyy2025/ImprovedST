# Computing Codes for the Paper "Improving the within-Node Estimation of Survival Trees while Retaining Interpretability"
### Haolin Li*, Yiyang Fan*, Jianwen Cai#


## Description

This repository contains computing codes for the paper "Improving the within-Node Estimation of Survival Trees while Retaining Interpretability". Please click [here](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.13821) for the full text of the paper.

## Folders

### 1-Data Generation

In this folder, we summarize the computing codes for generating survival data following a Weibull hazard model. The names of the code and the corresponding simulation scenarios in the paper are as follows,

## Continuous

# Low Dimension
* *data_generation_1.r* - tree-based effect, 100 sample size, 80% event rate.
* *data_generation_2.r* - tree-based effect, 100 sample size, 50% event rate.
* *data_generation_3.r* - tree-based effect, 300 sample size, 50% event rate.
* *data_generation_4.r* - tree-based effect, 300 sample size, 30% event rate.
* *data_generation_5.r* - linear covariate effect, 100 sample size, 80% event rate.
* *data_generation_6.r* - linear covariate effect, 100 sample size, 50% event rate.
* *data_generation_7.r* - linear covariate effect, 300 sample size, 50% event rate.
* *data_generation_8.r* - linear covariate effect, 300 sample size, 30% event rate.
* *data_generation_9.r* - nonlinear covariate effect, 100 sample size, 80% event rate.
* *data_generation_10.r* - nonlinear covariate effect, 100 sample size, 50% event rate.
* *data_generation_11.r* - nonlinear covariate effect, 300 sample size, 50% event rate.
* *data_generation_12.r* - nonlinear covariate effect, 300 sample size, 30% event rate.

### 2-Analysis

In this folder, we summarize the computing codes for fitting the proposed survival tree ensemble method and comparing the concordance index with single tree approach. The names and descriptions of the files are as follows,

* *analysis_1.r* - tree-based effect, 100 sample size, 80% event rate.
* *analysis_2.r* - tree-based effect, 100 sample size, 50% event rate.
* *analysis_3.r* - tree-based effect, 300 sample size, 50% event rate.
* *analysis_4.r* - tree-based effect, 300 sample size, 20% event rate.
* *analysis_5.r* - linear covariate effect, 100 sample size, 80% event rate.
* *analysis_6.r* - linear covariate effect, 100 sample size, 50% event rate.
* *analysis_7.r* - linear covariate effect, 300 sample size, 50% event rate.
* *analysis_8.r* - linear covariate effect, 300 sample size, 20% event rate.
* *analysis_9.r* - nonlinear covariate effect, 100 sample size, 80% event rate.
* *analysis_10.r* - nonlinear covariate effect, 100 sample size, 50% event rate.
* *analysis_11.r* - nonlinear covariate effect, 300 sample size, 50% event rate.
* *analysis_12.r* - nonlinear covariate effect, 300 sample size, 20% event rate.

### 3-Applications

In this folder, we include the datasets and computing code used to illustrate the application of the proposed method.

* *lung.csv* - a cleaned version of North Central Cancer Treatment Group Lung Cancer Data, which is publically available in the R package `survival' (Therneau and Lumley, 2013).
* *Full Data Tree Structure_lung.r* - the computing code for analyzing North Central Cancer Treatment Group Lung Cancer Data.
* *S1Data.csv* - a cleaned version of Cardiovascular Medical Records from the Faisalabad Institute of Cardiology, which is publically available (Chicco and Jurman, 2020).
* *Full Data Tree Structure_cardio.r* - the computing code for analyzing Cardiovascular Medical Records from the Faisalabad Institute of Cardiology.

## References

* Li, H., Fan, Y., & Cai, J. (2024+). Improving the within-Node Estimation of Survival Trees while Retaining Interpretability. Manuscript Submitted for Publication.
* Therneau, T., & Lumley, T. (2013). R survival package. R Core Team, 523.
* Chicco, D., & Jurman, G. (2020). Machine learning can predict survival of patients with heart failure from serum creatinine and ejection fraction alone. BMC medical informatics and decision making, 20, 1-16.
