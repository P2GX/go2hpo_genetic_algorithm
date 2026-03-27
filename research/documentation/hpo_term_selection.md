# HPO term candidates for Boolean GA

## Data snapshot
- Source files: `data/hpo2gene/phenotype_to_genes.txt`, GOA `data/gaf/goa_human.gaf.gz`, GTEx median TPM `data/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz`, HGNC map `data/map/hgnc_complete_set.txt`.
- Parsed 11,578 HPO terms with 5,128 unique gene symbols; 2,937 terms fall in the 20–200 gene window.
- GO and GTEx coverage is excellent: 5,055/5,128 (98.6%) HPO genes map to both datasets, so coverage is not a limiting factor.
- Baseline term HP:0001083 (Ectopia lentis): 46 genes (all covered by GO/GTEx).

## Filtering heuristic (practical search space)
- Focus on terms with 20–150 annotated genes to balance signal vs. over-specificity.
- Prefer phenotypes with recognizable, process-level GO signal (developmental pathways, organ-specific expression) over very broad metabolic panels.
- Keep a few very small (≈20 genes) terms where specificity may yield clean rules.

## Suggested candidate pools (by size band, all have full GO/GTEx coverage)
- Mid (≈60–80 genes, likely sweet spot):  
  HP:0004756 Ventricular tachycardia (70), HP:0002616 Aortic root aneurysm (71), HP:0001199 Triphalangeal thumb (70), HP:0010931 Abnormal blood sodium concentration (70), HP:0003215 Dicarboxylic aciduria (70), HP:0008665 Clitoral hypertrophy (70), HP:5200029 Social disinhibition (70).
- Small (≈20 genes, high specificity):  
  HP:0001022 Albinism (20), HP:0007754 Macular dystrophy (20), HP:0007350 Upper limb hyperreflexia (20), HP:0003387 Decreased number of large peripheral myelinated nerve fibers (20).
- Large but still ≤200 genes (for broader rules):  
  HP:0000664 Synophrys (199), HP:0009803 Short phalanx of finger (200), HP:0012622 Chronic kidney disease (197), HP:0000958 Dry skin (198).

## Quick GA probes (pop 40, 12 gens, mutation 0.45, β≈3 auto-estimated)
| HPO term | Genes | Best score (Fβ/penalty) | Best precision | Best recall | Notes |
| --- | --- | --- | --- | --- | --- |
| HP:0001022 Albinism | 20 | 0.59 | 0.52 | 0.60 | Melanosome organization term dominates; short rule. |
| HP:0001199 Triphalangeal thumb | 70 | 0.38 | 0.08 | 0.64 | Translation/RNA-binding + chromatin combo; recall-heavy. |
| HP:0004756 Ventricular tachycardia | 70 | 0.36 | 0.21 | 0.39 | Heart LV up-expression + Ca-channel terms. |
| HP:0001083 Ectopia lentis (baseline) | 46 | 0.28 | 0.13 | 0.33 | Remains modest; mixed ECM/developmental terms. |
| HP:0002616 Aortic root aneurysm | 71 | 0.24 | 0.08 | 0.31 | Weak signal; diffuse GO terms. |

Interpretation: small, specific phenotypes (albinism) yield substantially better Fβ; mid-size (≈70 genes) give moderate recall but weak precision. HP:0001083 remains challenging relative to these alternatives.

## Recommendations / next steps
- Prioritize runs on: HP:0001022, HP:0001199, HP:0004756, HP:0003215, HP:0010931. Treat HP:0002616/HP:0001083 as harder baselines.
- For tough terms (e.g., HP:0001083), try: lower β (1.5–2) to rebalance precision, add mild length penalty (`penalty_lambda` 0.05–0.1), and increase `max_n_terms` to 6–7 while keeping `max_n_conj` ≤4.
- For small terms (≈20 genes), keep population modest (20–40) but extend generations to 20–25 to refine precision without overgrowth; consider reducing mutation rate to 0.35.
- Capture stats to CSV via `--output-file` (already in `stats/` for the tested runs) and compare best-one precision/recall alongside rule size for convergence diagnostics.
- If precision remains low, consider constraining GO pools to phenotype-specific GO terms (`filtered_go_terms`) and enabling tissue priors (Heart/LV for arrhythmia, ocular tissues for lentis).

