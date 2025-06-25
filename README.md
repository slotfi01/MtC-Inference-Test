# Inference Test for Mean-to-CVaR (or Omega) financial performance ratios
Inference Test for Mean-to-CVaR (or Omega) financial performance ratios

This repository contains the code (in R and MATLAB) that implements the asymptotic inference test developed by 

Lotfi, S., Pagliardi, G., Paparoditis, E., & Zenios, S. A. (2025).  
Hedging political risk in international portfolios. European Journal of Operational Research.  
[https://doi.org/10.1016/j.ejor.2024.10.017](https://doi.org/10.1016/j.ejor.2024.10.017)

Bib citation: 

@article{LPZP2025,
title = {Hedging political risk in international portfolios},
journal = {European Journal of Operational Research},
volume = {322},
number = {2},
pages = {629-646},
year = {2025},
url = {https://www.sciencedirect.com/science/article/pii/S0377221724007938},
author = {Somayyeh Lotfi and Giovanni Pagliardi and Efstathios Paparoditis and Stavros A. Zenios}
}





## 📄 About the Paper

This paper shows how political risk —a global, non-diversifiable P-factor— affects international portfolio performance. It achieves four things:
- Develop a portfolio selection model that incorporates political risk hedging using a second-order stochastic dominance consistent criterion, the **Mean-to-CVaR (MtC) ratio.
- Introduce a novel asymptotic inference test for comparing MtC performance across portfolios, filling a methodological gap in the literature.
- Empirically show that politically hedged portfolios outperform both market indices and equally weighted portfolios across the U.S., Eurozone, and Japan.
- Demonstrate that political risk remains material even after currency hedging and during events like **COVID-19 and the war in Ukraine.






## 📂 What's in This Repository

- `mtc_cvr_inf_test.R` – R script to compute the MtC ratio and perform the inference test  
- `mtc_cvr_inf_test.m` – MATLAB script with equivalent implementation  
- `README.md` – Documentation of the project  
- `LICENSE` – License terms (MIT recommended)  
- `CITATION.cff` – Citation metadata for GitHub


## 🚀 Getting Started

mtc_cvr_inf_test test
                        H0: MTC(rt1) - MTC(rt0) = 0
                        H1: MTC(rt1) - MTC(rt0) > 0
 and
                        H0: CVaR(rt1) - CVaR(rt0) = 0
                        H1: CVaR(rt1) - CVaR(rt0) > 0
 and
                        H0: SHR(rt1) - SHR(rt0) = 0
                        H1: SHR(rt1) - SHR(rt0) > 0

and report the p-values along with differences of MTCs and CVaRs

Inputs:
b: Length of block
B: Number of bootstrap iterations
rt1, rt0: Return (or excess return) time series with length of S from
which we resample blocks of length b with replacement
alpha: Confidence level for calculating CVaR and MTC

Outputs:
mtc_diff: MTC(rt1) - MTC(rt0)
cvr_diff: CVaR(rt1) - CVaR(rt0)
mtc_pval: The p-value of testing MTC differences
cvr_pval: The p-value of testing CVaR differences




