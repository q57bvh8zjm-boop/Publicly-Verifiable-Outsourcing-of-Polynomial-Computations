# Benchmark configuration

Benchmark records use the following machine-readable CSV fields:

```text
scheme,d,t,m,trial,Tn_sec,Tp_sec,Tv_sec,Td_sec,Ts_sec,Rt_star,log2_Rt_star
```

Measurements use a monotonic wall-clock timer. KeyGen is recorded separately and excluded from $Rt^*$, matching the paper definition.

## Parameter grids

- Figure 2: $d=2,t=1,m=200,400,\ldots,2000$ for $\Gamma_1,\Gamma_2,\Gamma_3$.
- Figure 3: $d=2,t=2,3,m=200,400,\ldots,2000$ for $\Gamma_1,\Gamma_2$.
- Figure 4: $d=4,8,10,t=1,m=1,\ldots,15$ for $\Gamma_1,\Gamma_2,\Gamma_3$.
- Figure 5: $d=2,3,t=1,2,m=200,400,\ldots,2000$ for $\Gamma_1$ and $\Pi_3$.
- Microbenchmark: 100 trials for $\mathsf{Add}_q$, $\mathsf{Mul}_q$, $\mathsf{Add}_{q^2}$, $\mathsf{Mul}_{q^2}$, $\mathsf{Exp}_{\mathbb G}$, and Paillier/HSS encryption.

Polynomial inputs use the coefficient order documented in the top-level README and allocate `mpoly_num_coefficients(m, d)` coefficients.

The corresponding result datasets are stored in `../results/reported/`.
