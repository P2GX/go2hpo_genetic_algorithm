 # go2hpo_genetic_algorithm
 
 Genetic Algorithm for discovering DNF logic formulas that explain HPO term - gene associations.
 
 ## Prerequisites
 - Rust toolchain (cargo)
 - Access to the bundled data files under `data/` (used by the binaries)
 
 ## Running the algorithm
 Three entrypoints are provided, depending on whether you want a single run, an interactive session, or a batch of runs from a config file.
 
 ### Single run: `ga`
 Execute one GA run with CLI flags.
 
```sh
cargo run --release --bin ga -- --hpo-term HP:0001083 -p 20 -g 5 -m 0.5 -t 3 --max-n-terms 5 --max-n-conj 4 --rng-seed 42
```
 
 Key options (all long flags unless noted):
 - `--hpo-term` (required) HPO term ID, e.g. `HP:0001083`
 - `-p, --pop-size` default 20
 - `-g, --generations` default 5
 - `-m, --mutation-rate` default 0.5
 - `-t, --tournament-size` default 3
 - `--max-n-terms` default 5 (max terms per conjunction)
 - `--max-n-conj` default 4 (max conjunctions per DNF)
 - `-l, --penalty-lambda` default 0.0
 - `-b, --fscore-beta` optional (if omitted, estimated from class imbalance)
 - `-o, --output-file` optional (basename, saved under `stats/`)
 - `--rng-seed` default 42
 - `--export-bin / --import-bin` optional bincode snapshot handling
 - `--use-expanded` use pre-expanded GO annotations
 - `--disallow-go-negations` forbid NOT(GO:...) in formulas
 - Enrichment controls (default constants in code):
   - `--use-enriched-go-pool`
   - `--go-enrichment-top-k`
   - `--go-enrichment-p-value`
   - `--go-enrichment-min-support`
   - `--go-enrichment-include-parents`
   - `--go-enrichment-min-fold`
   - `--go-enrichment-max-bg-freq`
   - `--go-enrichment-filter-roots`
 
 ### Interactive session: `ga_interactive`
 Guided prompts let you run multiple HPO terms in one session.
 
```sh
cargo run --release --bin ga_interactive
```
 
 Flow:
 1) Choose GO annotation mode (expanded vs direct), GO negation allowance, and optional GO enrichment plus its thresholds.
 2) Loop: enter an HPO term, then per-run settings (population, generations, mutation rate, tournament size, max terms/conj, penalty, optional output file).
 3) Beta is estimated from class imbalance; you can accept or override.
 4) Type `quit` to exit.
 
 ### Batch runs from JSON: `ga_batch`
 Run one or more GA configs from a JSON list.
 
```sh
cargo run --release --bin ga_batch -- src/bin/run_configs/hpo_focused_one_run.ljson \
  --summary-file hpo_focused_one_run.txt \
  --summary-csv stats/runs/summaries/hpo_focused_one_run.csv
```
 
 Arguments:
 - `<config_path>` (positional, optional) JSON array of run objects; defaults to `src/bin/run_configs/hpo_runs.ljson`.
- `--summary-file, -s <path>` optional text table (relative paths are automatically placed under `stats/runs/summaries/`).
 - `--summary-csv, -c <path>` optional CSV of best-per-run stats.
 
 Notes:
 - Each run object maps to `GaConfig` (see `src/bin/ga_batch.rs` for fields). Batch mode always uses pre-expanded GO annotations.
 - Example config `src/bin/run_configs/hpo_focused_one_run.ljson` contains a single run for `HP:0005938`.
- All example commands above use release mode via `cargo run --release ...`.
 
