# GA runtime/parallelization notes

## Hotspots observed in code
- Fitness evaluation: every call to `ConjunctionScorer::fitness` / `DNFScorer::_fitness` rebuilds a satisfaction map by scanning all genes (`checker.all_satisfactions`). In `DNFScorer`, this happens once per conjunction and merges via OR; all serial and HashMap-heavy.
- Satisfaction checking: `NaiveSatisfactionChecker::are_term_annotations_satisfied` walks GO descendants per term for every gene check; `all_satisfactions` iterates all genes for each conjunction, sequentially.
- GA loop: `evolve_one_generation` builds offspring serially (selection → crossover → mutation → evaluation). Initial population generation (`new_with_size`) is also serial.
- Random GO term sampling uses `iter_term_ids().nth(rnd_index)` (O(n) walk) for every add/mutate.
- Logging: per-generation `println!` in `fit_with_stats_history` adds overhead when many generations/runs are batched.

## Parallelization opportunities
- Satisfaction per gene: `all_satisfactions` can use `rayon::prelude::*` to `par_iter()` the gene map (genes independent). Return could stay as `HashMap` or be switched to `Vec<bool>` keyed by stable gene indices.
- DNF satisfaction merge: represent conjunction satisfactions as bitsets (bitvec) and parallel-OR across conjunctions; avoids HashMap updates and enables SIMD-ish OR.
- Population evaluation: after generating offspring formulas, evaluate fitness in parallel (`par_iter_mut`) since `FormulaEvaluator` is immutable; similarly, initial population generation can `par_iter` over indices.
- Reproduction: selection + mutation are RNG-heavy; keep them serial or give each worker its own `SmallRng` seeded from a master seed to safely parallelize parts of offspring creation if needed.

## Quick wins (low risk)
- Precompute GO term pool into `Vec<TermId>` once; replace `iter_term_ids().nth` with `choose` on the vec (O(1) indexing).
- Add a verbosity flag to suppress per-generation logging during sweeps.
- Reuse buffers: avoid reallocating populations by swapping vectors; avoid cloning formulas when not needed (e.g., store score separately and mutate in place).
- Cache GO ancestor/descendant closures (e.g., `HashMap<TermId, SmallVec<TermId>>`) to remove repeated `iter_descendant_ids` walks in satisfaction checks.

## Higher-impact changes (moderate effort)
- Replace `HashMap<GeneId, bool>` with dense `Vec<bool>` keyed by an index map (`GeneId` → usize) stored in `GeneSetAnnotations`; this speeds both `all_satisfactions` and DNF OR merges.
- Add a parallel `all_satisfactions_par` and wire it into `ConjunctionScorer`/`DNFScorer`; keep the naive version for compatibility.
- Parallelize evaluation in `evolve_one_generation`: generate offspring formulas (serial RNG), then `par_iter_mut` to score them; likewise for `initialize_population` / `new_with_size`.
- Optional: parallelize per-generation stats (min/avg/max) with rayon reductions—small gain but easy once `rayon` is in place.

## RNG handling for parallel sections
- Use thread-local `SmallRng` seeds derived from a master seed (`seed ^ thread_idx as u64`) when parallelizing any mutation/selection path.
- Keep evaluator/scorer deterministic given the seed (no RNG inside fitness), enabling safe parallel scoring.

## Suggested implementation order
1) Add rayon dependency; precompute GO term pool vecs; guard logging with a verbosity flag.
2) Introduce dense gene index map + bitset satisfaction; parallelize `all_satisfactions`.
3) Parallelize fitness evaluation in population init and per-generation scoring; keep reproduction RNG serial at first.
4) If needed, parallelize reproduction with per-thread RNGs; benchmark vs. serial to ensure gains.

