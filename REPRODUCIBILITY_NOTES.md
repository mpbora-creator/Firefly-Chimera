# Reproducibility Notes

This package is the current-code version for Mac/Linux validation.

The exact validation case is:

```text
P = 9
kappa = 0.0485
tau0 = 1.185
sigtau = 0.062
ell = 12.0
delay_seed = 29
init_seed = 202
```

Run:

```bash
julia --project=. scripts/validate_exact_case.jl
```

Expected current-code result:

```text
score = 4.8783576387658645
chorus_fraction = 0.320625
chorus_count = 513
char_count = 1087
freq_mean = 0.38068713011729466
freq_std = 0.0045542274862573765
freq_cv = 0.011963176913425362
limit_cycle_center = (0.3950931921372677, 7.536484564072712)
limit_cycle_points = 6000
```

Important: resume from an existing live CSV is disabled by default.  This
prevents older classifications from being silently reused.  To intentionally
resume from a live CSV, run with:

```bash
CHIMERA_RESUME=1 julia --project=. scripts/search.jl
```

For cross-machine validation, use a clean package directory and do not copy old
`results/` files from previous diagnostic versions.
