# Chimera states in a network of firefly oscillators

This repository contains the Julia simulation, analysis, validation, and
figure-generation code accompanying the manuscript *Chimera states in a
network of firefly oscillators* by Deeptarag Das, Debashis Saikia, Kalpana
Bora, and Madhurjya P. Bora.

The code studies identical nonlocally coupled relaxation oscillators with
uniform or heterogeneous communication delays. The deposited workflows focus
on the heterogeneous-delay Sarfati–Peleg-like state, its chorus/character
diagnostics, and the reduced analytical description.

## System requirements

- Julia 1.12 (the deposited `Manifest.toml` was generated with Julia 1.12.6)
- A macOS, Linux, or Windows computer
- Sufficient memory and runtime for the requested simulation; the exact-case
  validation and parameter searches are substantially more expensive than
  package setup

No proprietary software is required. Julia dependencies are fixed by
`Manifest.toml`.

## Installation

Clone or download the repository, enter its root directory, and instantiate
the environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Precompile once if desired:

```bash
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

## Reproduce the reference case

Run the deterministic numerical check from the repository root:

```bash
julia --project=. scripts/validate_exact_case.jl
```

Success is reported as `VALIDATION PASSED`. The parameter set, expected
diagnostics, and cross-machine cautions are recorded in
[`REPRODUCIBILITY_NOTES.md`](REPRODUCIBILITY_NOTES.md).

## Main workflows

Generate diagnostics from completed paper cases:

```bash
CHIMERA_MODE=paper julia --project=. -t auto scripts/paper.jl
```

Search the parameter space and then run robustness calculations:

```bash
julia --project=. -t auto scripts/search.jl
CHIMERA_MODE=robustness julia --project=. -t auto scripts/robustness.jl
```

Search mode writes `results/natural_chimera_best_so_far.txt`. Paper mode reads
`results/natural_chimera_results_live.csv`. Input locations can be overridden:

```bash
ROBUST_BEST_FILE=my_best.txt CHIMERA_MODE=robustness julia --project=. scripts/robustness.jl
PAPER_SOURCE_CSV=my_results.csv CHIMERA_MODE=paper julia --project=. scripts/paper.jl
```

Selected reduced-model validation workflows are:

```bash
CHIMERA_QUICK=1 julia --project=. scripts/model_validation.jl
CHIMERA_QUICK=1 julia --project=. scripts/analytical_validation.jl
CHIMERA_QUICK=1 julia --project=. scripts/revised_analytical_validation.jl
CHIMERA_QUICK=1 julia --project=. scripts/distributed_delay_validation.jl
CHIMERA_QUICK=1 julia --project=. scripts/heterogeneous_state_validation.jl
CHIMERA_QUICK=1 julia --project=. scripts/argument_principle_validation.jl
CHIMERA_QUICK=1 julia --project=. scripts/reduction_robustness_validation.jl
CHIMERA_QUICK=1 julia --project=. scripts/time_dependent_reduced_validation.jl
```

`CHIMERA_QUICK=1` performs a reduced validation intended for checking the
workflow. Omit it for the full research calculation.

## Repository contents

- `src/`: simulation, diagnostics, plotting, and reduced-model modules
- `scripts/`: executable search, validation, and figure entry points
- `results/`: deposited numerical results and panel-level plotting data
- `figures/`: validation figures and their underlying tabular data
- `Project.toml`, `Manifest.toml`: Julia environment specification
- `FIGURE_CODE_MAP.md`: mapping from manuscript results to commands and files
- `CODE_AVAILABILITY.md`: manuscript statement and release placeholders
- `CITATION.cff`: machine-readable software citation metadata

See [`FIGURE_CODE_MAP.md`](FIGURE_CODE_MAP.md) before attempting individual
manuscript panels. The deposited panel-level `.dat` and `.csv` files allow the
reported plots to be inspected without rerunning long searches.

## Reproducibility notes

Random seeds for the reference case are fixed. Resume from an existing live
CSV is disabled by default so that results from older diagnostic definitions
are not silently reused. To resume intentionally:

```bash
CHIMERA_RESUME=1 julia --project=. scripts/search.jl
```

For a clean cross-machine check, do not copy old `results/` checkpoint files
into a fresh clone. Small floating-point differences may occur across Julia
versions or platforms; use Julia 1.12 and the deposited manifest for the
closest reproduction.

## Citation and license

Before publication, cite the archived Zenodo release using the metadata in
`CITATION.cff`. Release instructions and manuscript wording are in
[`CODE_AVAILABILITY.md`](CODE_AVAILABILITY.md).

The source code is released under the [MIT License](LICENSE). Deposited article
text, figures, and third-party material may be subject to their own publication
terms.
# Firefly-Chimera
