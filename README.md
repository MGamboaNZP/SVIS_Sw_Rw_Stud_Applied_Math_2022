# Numerical computation for moments and distribution of Sw and Rw in a stochastic SVI model

This repository contains the MATLAB implementation used to compute:

- The **distribution of \(S_w\)**  
- The **distribution of \(R_w\)**  
- The **first and second moments** associated with these random variables  

in the **stochastic SVI (Susceptible–Vaccinated–Infected)** epidemic model with imperfect vaccine, as introduced in:

> Gamboa, M., & Lopez‐Herrero, M. (2022).  
> *Measures to assess a warning vaccination level in a stochastic SIV model with imperfect vaccine.*  
> Studies in Applied Mathematics, 148(4), 1411–1438.  
> https://doi.org/10.1111/sapm.12479

## Zenodo doi
https://doi.org/10.5281/zenodo.18093320

The algorithms implemented here reproduce the numerical procedures described in the paper and extend them with **stable, vectorized and parallelized** routines for large population sizes.

---

## 🚀 Features

- Computation of the **distribution of \(S_w\)** for arbitrary warning levels \(w\)
- Computation of the **distribution of \(R_w\)** (if included in your repo)
- Stable numerical implementation using:
  - log‑space arithmetic  
  - vectorized operations  
  - optional parallelization (`parfor`)
- Reproducible MATLAB scripts corresponding to the published results
- Ready-to-use code for large-scale simulations (e.g., \(N = 500\), \(v_0 = 480\))

---

## 📂 Repository structure

