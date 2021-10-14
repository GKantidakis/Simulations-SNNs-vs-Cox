# Simulations-SNNs-vs-Cox

This repository stores the R-code of a **simulation study** to compare survival neural networks (SNNs) with Cox models for clinical trial data.  

The predictive performance of ML techniques is compared with statistical models in a simple clinical setting (small/moderate sample size, small number of predictors) with Monte Carlo simulations. Synthetic data (250 or 1000 patients) are generated that closely resemble 5 prognostic factors pre-selected based on a European Osteosarcoma Intergroup study (MRC BO06/EORTC 80931).  

Comparison is performed between two partial logistic artificial neural networks (PLANN original by Biganzoli et al. 1998, Statistics in medicine, 17(10), 1169-1186 and PLANN extended by Kantidakis et al. 2020 BMC medical research methodology, 20(1), 1-14) as well as Cox models for 20, 40, 61, and 80% censoring. Survival times are generated from a log-normal distribution. Models are contrasted in terms of C-index, Brier score at 0-5 years, Integrated Brier Score (IBS) at 5 years, and miscalibration at 2 and 5 years. Endpoint of interest is overall survival.  

Note: PLANN original/extended are tuned based on IBS at 5 years and C-index.
