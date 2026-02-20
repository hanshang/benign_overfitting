# Code Repository for Benign Overfitting in High-Dimensional Linear Discriminant Analysis

## Folder structure

```
kehan_paper/
├── README.md
├── code/
│   ├── run_all.R              # Master script that runs everything and generates all figures
│   ├── ARCENE/                # Raw ARCENE cancer dataset (train/valid/test splits and labels)
│   ├── inputs/                # Generated CSV files from simulations and ARCENE processing
│   └── scripts/               # R scripts for the analysis pipeline
│       ├── experiments_lib.R        # Shared helpers (LDA fitting, sampling, covariance constructors, theory formulas)
│       ├── generate_inputs.R        # Monte Carlo simulations and analytic theory curves
│       ├── generate_arcene_inputs.R # ARCENE dataset processing and classifier comparisons
│       ├── generate_figures.R       # Reads CSVs from inputs/ and writes PNGs to plots/
│       └── generate_missing_inputs.R# Re-generates any missing input CSVs
└── plots/                     # Output figures (PNGs) used in the paper
```

## One-command run

```sh
Rscript code/run_all.R --clean --paper
```

`--paper` uses 100 Monte Carlo repetitions per curve point.

Other flags:

- `--clean` -- delete all generated inputs and plots before running
- `--force` -- recompute all CSVs even if they already exist
- `--quick` -- use 8 reps (fast, for testing)
- `--seed=N` -- override the default random seed (123)
- `--reps=N` -- set a custom number of MC repetitions

## Dependencies

- Base R + `MASS`
