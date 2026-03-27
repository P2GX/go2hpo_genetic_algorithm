## GO Enrichment Pool Flags

This feature constrains GO-term sampling to a per-HPO enriched pool, instead of the full ontology or the simple union of all GO terms annotated to positive genes. It reduces noise for phenotypes with many positives and speeds convergence toward higher-precision rules.

- **`--use-enriched-go-pool`** (bool, default `false`)  
  When true, build a phenotype-specific GO pool via enrichment before running the GA. If the enriched pool is empty, the code falls back to the union of positive genes’ GO terms; if that is empty, it falls back to the full GO ontology.

- **`--go-enrichment-top-k`** (usize, default `200`)  
  Keep at most this many enriched GO terms after sorting by p-value (and fold-enrichment tie-break). Limits search-space size and runtime.

- **`--go-enrichment-p-value`** (f64, default `0.05`)  
  Right-tailed hypergeometric (Fisher) p-value cutoff. Only terms with `p <= cutoff` are retained before the top-k filter. Lower values tighten the pool; higher values broaden it.

- **`--go-enrichment-min-support`** (usize, default `2`)  
  Minimum number of positive genes annotated with a GO term to consider it for enrichment. Filters out very rare/unstable terms.

- **`--go-enrichment-include-parents`** (bool, default `true`)  
  If enabled, add direct parents of enriched terms to the pool to keep some generalization and avoid over-fragmenting on very specific child terms.

- **`--go-enrichment-min-fold`** (f64, default `1.2`)  
  Keep terms only if fold >= min_fold or fold <= 1/min_fold. Symmetric, so it retains enriched or depleted terms and drops near-1.0 roots/ubiquitous terms.

- **`--go-enrichment-max-bg-freq`** (f64, default `0.2`)  
  Skip GO terms whose background frequency exceeds this fraction of all genes. Set to `1.0` to disable. Removes overly common/low-information terms.

- **`--go-enrichment-filter-roots`** (bool, default `true`)  
  Drop GO root/namespace-level terms (e.g., GO:0008150, GO:0003674, GO:0005575, GO:0009987, GO:0110165).

### Why these flags exist
- **Noise reduction:** The prior pool (union of positives) can be large and noisy for common phenotypes. Enrichment trims low-signal terms.
- **Faster convergence:** Smaller, phenotype-relevant search space improves GA efficiency and can raise best F-scores.
- **Safety rails:** Top-k, p-value, min-support, min-fold, root/high-prevalence filters, and parent inclusion keep the pool informative without exploding or vanishing.
- **Backwards compatibility:** Flags default to the previous behavior (enrichment off), so existing configs keep working.

### Example command
Run the standard GA binary with enrichment on (using defaults for other enrichment knobs):
```
cargo run --release --bin ga -- --hpo-term HP:0001083 --use-expanded --use-enriched-go-pool -p 140 -g 50
```

