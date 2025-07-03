# GBS-Enhanced Initialization for Jaccard-Weighted Densest Subgraph Problem

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Status](https://img.shields.io/badge/status-exploratory-lightgrey)

> **Author:** Yashaswini Mathur  
> **Affiliations:**  
> - TCS Research & Innovation, Bangalore  
> - IISER Bhopal  

---

## Overview

This project explores a **quantum-enhanced hybrid approach** to solving the **Jaccard-Weighted Densest Subgraph (JWDS)** problem, a challenging NP-hard extension of the classical Densest Subgraph Problem (DSP).

We integrate **Gaussian Boson Sampling (GBS)**, a photonic quantum sampling method, into the **initialization phase** of the classical **Iterative (ITR)** algorithm to assess potential improvements in optimization performance.

The codebase contains:
- Classical JWDS-ITR algorithm implementations
- GBS-based initialization subroutines
- Evaluation metrics: density, Jaccard similarity, and total Q-score
- Example graphs and test runs

This work was conducted as part of a research internship at **TCS Research**, aiming to demonstrate how **quantum technologies** can offer practical advantages for real-world classical approximation problems.

---

## Report

*[Internship Project Report (PDF)](link-to-report-if-hosted)*

> **Title**: *Exploration of Quantum-Enhanced Initialization for Jaccard-Weighted Densest Subgraphs*  
> This report presents a detailed study of the JWDS problem, the classical ITR algorithm, and a proposed GBS-assisted hybrid approach. It outlines:
> - Theoretical foundations of DSP and JWDS  
> - Overview of the ITR algorithm  
> - Introduction to Gaussian Boson Sampling  
> - Classical and GBS-based implementations  
> - Preliminary performance evaluations  
> - Future simulation plans

---

## Codes
This repo consists of the following code files:
1. `Weighted_Jaccard_Iterative_Algo2.py` -  This the benchmark ITR algorithm from Arachchi et al. (2024)
2. `Weighted_Jaccard_ITR_with_GBS.py` - This is the enhanced GBS variation we introduced
3. `basics_of_gbs(pennylane).ipynb` - This is the code for Appendix A
4. `dense_subgraph.ipynb` - This is the code for Appendix A
5. `training_variational_gbs_distribuitions.ipynb` - This is the code for Appendix B
6. `Jaccard_Iterative_Algo1.py` - This is the initial implementation of Algo 1 from Arachchi et al. (2024)
7. `Weighted_Jaccard_Iterative_Algo3.py` - This the greedy algorithm from Arachchi et al. (2024) 

**Note:** To run this code, install the dependencies:
```bash
pip install networkx numpy strawberryfields
```

---

## Other files
1. `AppendixA_Introduction_and_Basics` - A walkthrough
2. `AppendixB_Training_Variational_GBS_Distributions` - A walkthrough
3. `AppendixC_Graph_Similarity` - A walkthrough

---

## Key Concepts
1. JWDS (Jaccard-Weighted Densest Subgraph): Extension of DSP where we maximize both intra-subgraph density and similarity across multiple graphs.
2. ITR Algorithm: A 2-approximation iterative algorithm that greedily refines vertex sets for each graph to improve the joint objective.
3. GBS (Gaussian Boson Sampling): A photonic quantum process shown to preferentially sample dense subgraphs — used here to generate high-quality initial vertex sets.

---

## Metrics Tracked
1. Density of each subgraph
2. Jaccard Index between subgraphs (pairwise similarity)
3. Q-Score: Overall objective combining density and Jaccard similarity, weighted by a tunable λ-parameter

---

## Future Work
This project sets the groundwork for:
1. Full GBS circuit simulation using Strawberry Fields
2. Comparing GBS-based initialization to classical baselines across various graph sizes
3. Robustness and convergence speed analysis
4. Studying other quantum-inspired or hybrid heuristics

