# Parallelizing the GA: safest path first

Context: `ga_common::run_ga` builds a `GeneticAlgorithm`, whose `evolve_one_generation` drives selection → crossover → mutation → fitness via `FormulaEvaluator` (see `src/genetic_algorithm/algorithm.rs`). Fitness calls flow into `DNFScorer` → `ConjunctionScorer` → `NaiveSatisfactionChecker::all_satisfactions` (see `src/logical_formula/satisfaction_checker.rs`), which scans every gene and walks GO descendants for each literal.

## Where the time goes
- Satisfaction is the hotspot: `all_satisfactions` is `O(#genes × #literals)` and walks GO descendants on every check; `DNFScorer` then ORs these per-conjunction results.
- GA loop is serial: offspring are generated and evaluated one-by-one inside `evolve_one_generation`, even though fitness has no RNG and is pure.
- Logging and iterator choices add smaller overhead (per-generation `println!`, repeated `iter_descendant_ids`, and random GO sampling via `iter_term_ids().nth`), but they are secondary to the above.

## Top vs bottom parallelization
- Bottom (inside fitness) is safer to start with: `NaiveSatisfactionChecker::all_satisfactions` has no shared mutability or RNG. Turning its gene iteration into `par_iter()` gives parallelism where the CPU time is spent, without touching GA control flow.
- Top (population-level) parallelism is possible but needs extra trait bounds (`FitnessScorer: Send + Sync`, `FormulaEvaluator` shared safely) and must respect RNG ownership. It also risks nested Rayon if fitness itself becomes parallel.
- Avoid nested Rayon: pick one layer. If fitness uses Rayon internally, keep population evaluation serial (or gate it behind a feature); if population evaluation goes parallel, keep `all_satisfactions` serial or move it to a non-Rayon bitset path.

## Recommended low-risk sequence
✅ 1) Add Rayon dependency and guard logging with a verbosity flag (free win for sweeps). Precompute GO term pools into `Vec<TermId>` to replace `iter_term_ids().nth`.
✅ 2) Parallelize satisfaction: in `NaiveSatisfactionChecker::all_satisfactions`, switch to `par_iter()` over `get_gene_annotations_map()` and collect into a `HashMap`. This keeps the public API stable and avoids RNG concerns. If/when needed, add a cached GO descendant map to cut repeated `iter_descendant_ids` calls.
3) Optional, after benchmarking step 2: evaluate offspring in parallel. Restructure `evolve_one_generation` into two phases—(a) build offspring formulas serially (RNG-safe); (b) score them with `par_iter_mut`, using a `&FormulaEvaluator` that is `Send + Sync`. Leave selection/mutation serial unless you introduce per-thread RNGs.
4) Later, consider a dense satisfaction mask (e.g., `bitvec` or `Vec<u8>`) keyed by a `GeneId → usize` map to speed up DNF OR merges. Avoid `Vec<bool>` because of its bit-packed accessor costs; prefer `bitvec` or byte masks for predictable access.

## Pitfalls to avoid
- No mutexes inside fitness; keep data read-only and precompute caches instead.
- Avoid branch-heavy inner loops (structure exclusion/inclusion checks so the hot path is straight-line).
- Do not nest Rayon (`par_iter` inside another `par_iter`) without bounding the thread pool.
- Skip `Vec<bool>` for hot paths; use `bitvec` or `Vec<u8>` if you move to dense masks.

