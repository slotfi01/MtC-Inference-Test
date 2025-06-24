# MtC-Inference-Test
inference test for mean-to-CVaR ratios

# Inference Test for Politically Hedged Portfolios

This repository contains the code (in R and MATLAB) that implements the asymptotic inference test developed in the paper:

Lotfinoghabi, S., Pagliardi, G., Paparoditis, E., & Zenios, S. A. (2025).  
Hedging political risk in international portfolios. European Journal of Operational Research.  
[https://doi.org/10.1016/j.ejor.2024.10.017](https://doi.org/10.1016/j.ejor.2024.10.017)


@article{LPZP2025,
title = {Hedging political risk in international portfolios},
journal = {European Journal of Operational Research},
volume = {322},
number = {2},
pages = {629-646},
year = {2025},
issn = {0377-2217},
doi = {https://doi.org/10.1016/j.ejor.2024.10.017},
url = {https://www.sciencedirect.com/science/article/pii/S0377221724007938},
author = {Somayyeh Lotfi and Giovanni Pagliardi and Efstathios Paparoditis and Stavros A. Zenios},
keywords = {Block bootstrapping, Conditional value-at-risk, Equity home bias puzzle, International diversification, Portfolio selection, Skwewed returns},
abstract = {We show that internationally diversified portfolios carry sizeable political risk premia and expose investors to tail risk. We obtain political efficient frontiers with and without hedging political risk using a portfolio selection model for skewed distributions and develop a new asymptotic inference test to compare portfolio performance. Politically hedged portfolios outperform a broad market index and the equally weighted portfolio for US, Eurozone, and Japanese investors. Political risk hedging is not subsumed by currency hedging, and the diversification gains of politically hedged portfolios persist under currency hedging and transaction cost frictions. Hedging political risk induces equity home bias but does not fully explain the puzzle.}
}



## 📄 About the Paper

This paper addresses how **political risk**—a global, non-diversifiable factor—affects international portfolio performance. The authors:

- Develop a portfolio selection model that incorporates political risk hedging using a second-order stochastic dominance (SSD) consistent criterion, the **Mean-to-CVaR (MtC) ratio.
- Introduce a novel asymptotic inference test for comparing MtC performance across portfolios, filling a methodological gap in the literature.
- Empirically show that politically hedged portfolios outperform both market indices and equally weighted portfolios across the U.S., Eurozone, and Japan.
- Demonstrate that political risk remains material even after currency hedging and during events like **COVID-19 and the war in Ukraine.





## 📂 What's in This Repository

- `inference_test.R` – R script to compute the MtC ratio and perform the inference test  
- `inference_test.m` – MATLAB script with equivalent implementation  
- `README.md` – Documentation of the project  
- `LICENSE` – License terms (MIT recommended)  
- `CITATION.cff` – Citation metadata for GitHub



