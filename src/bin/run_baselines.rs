use go2hpo_genetic_algorithm::{
    annotations::GeneSetAnnotations,
    genetic_algorithm::ScoreMetric,
    other_models::{
        run_constant_baseline, ConstantStrategy,
        best_single_term, SingleTermConfig,
        greedy_forward_rule, GreedyConfig,
    },
    logical_formula::PreexpandedSatisfactionChecker,
    utils::fixtures::gene_set_annotations::{gene_set_annotations_expanded, go_ontology},
};
use ontolius::TermId;
use std::env;

fn main() -> anyhow::Result<()> {
    // 1) Load ontology and an expanded gene set (uses cache if present)
    let go = go_ontology();
    let gene_set: GeneSetAnnotations = gene_set_annotations_expanded(&go);

    // 2) Choose the phenotype to evaluate (arg0 = binary, arg1 = HPO term). Default: HP:0000001
    let phenotype: TermId = env::args()
        .nth(1)
        .unwrap_or_else(|| "HP:0000001".to_string())
        .parse()
        .expect("Invalid HPO term ID");

    // 3) Choose a recall-heavy metric
    let metric = ScoreMetric::FScore(3.0);

    // 4) Build a satisfaction checker over the expanded annotations (fast, no ontology lookups)
    let checker = PreexpandedSatisfactionChecker::new(&gene_set);

    // 4) Constant baseline (all-positive vs all-negative; pick best by metric)
    let const_res =
        run_constant_baseline(&gene_set, &phenotype, &metric, ConstantStrategy::MetricBest);
    println!("Constant: {:?}", const_res);

    // 5) Single-term baseline (best literal)
    let st_conf = SingleTermConfig {
        metric: metric.clone(),
        min_support: 1,           // bump if you want to ignore ultra-rare hits
        min_improvement: 1e-6,    // metric tolerance
        min_recall_gain: 0.01,    // require recall gain on ties
    };
    let st_res = best_single_term(&checker, &phenotype, &st_conf);
    println!("Single-term: {:?}", st_res);

    // 6) Greedy forward rule-builder
    let g_conf = GreedyConfig {
        metric,
        max_literals: 5,          // cap size
        min_improvement: 1e-6,    // metric tolerance
        min_support: 1,
        min_recall_gain: 0.01,    // allow recall-driven steps
    };
    let g_res = greedy_forward_rule(&checker, &phenotype, &g_conf);
    println!("Greedy: {:?}", g_res);

    Ok(())
}