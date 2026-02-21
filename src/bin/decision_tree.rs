use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use go2hpo_genetic_algorithm::{
    annotations::GeneSetAnnotations,
    genetic_algorithm::ScoreMetric,
    logical_formula::{NaiveSatisfactionChecker, PreexpandedSatisfactionChecker},
    other_models::{decision_tree_rules, DecisionTreeConfig, DecisionTreeLeaf, DecisionTreeResult},
    utils::fixtures::gene_set_annotations::{
        gene_set_annotations, gene_set_annotations_expanded, go_ontology,
    },
};
use ontolius::ontology::{csr::MinimalCsrOntology, HierarchyWalks, OntologyTerms};
use ontolius::term::MinimalTerm;
use ontolius::TermId;
use statrs::distribution::{DiscreteCDF, Hypergeometric};

// Defaults aligned with GA runner.
const DEFAULT_GO_ENRICHMENT_TOP_K: usize = 200;
const DEFAULT_GO_ENRICHMENT_MIN_SUPPORT: usize = 2;
const DEFAULT_GO_ENRICHMENT_P_VALUE: f64 = 0.05;
const DEFAULT_GO_ENRICHMENT_INCLUDE_PARENTS: bool = true;
const DEFAULT_GO_ENRICHMENT_MIN_FOLD: f64 = 1.2;
const DEFAULT_GO_ENRICHMENT_MAX_BG_FREQ: f64 = 0.2;
const DEFAULT_GO_ENRICHMENT_FILTER_ROOTS: bool = true;

/// Train a shallow decision tree baseline and print the learned conjunction rules (DNF).
#[derive(Parser, Debug)]
#[command(name = "decision_tree")]
#[command(about = "Decision-tree baseline that outputs conjunction rules for an HPO term")]
struct Args {
    /// HPO term to analyze (e.g., HP:0001083)
    #[arg(long)]
    hpo_term: String,

    /// Maximum tree depth (root = 0)
    #[arg(long, default_value_t = 4)]
    max_depth: usize,

    /// Minimum number of genes a rule must cover to be considered
    #[arg(long, default_value_t = 1)]
    min_support: usize,

    /// Minimum improvement needed to accept a split
    #[arg(long, default_value_t = 1e-6)]
    min_improvement: f64,

    /// Minimum score for a leaf to be kept as a rule
    #[arg(long, default_value_t = 1e-6)]
    min_leaf_score: f64,

    /// Metric: one of accuracy|precision|recall|f|f1|f2|f0.5|f{beta}
    #[arg(long, default_value = "f1")]
    metric: String,

    /// F-score beta (used when metric starts with 'f')
    #[arg(long, default_value_t = 1.0)]
    beta: f64,

    /// Use pre-expanded GO annotations (direct + ancestors)
    #[arg(long, default_value_t = false)]
    use_expanded: bool,

    /// Disallow NOT(GO:...) literals in rules
    #[arg(long, default_value_t = false)]
    disallow_go_negations: bool,

    /// Use enriched GO pool (hypergeometric top-k) instead of simple union of positives
    #[arg(long, default_value_t = false)]
    use_enriched_go_pool: bool,

    /// Max number of enriched GO terms to keep (before optional parent expansion)
    #[arg(long, default_value_t = DEFAULT_GO_ENRICHMENT_TOP_K)]
    go_enrichment_top_k: usize,

    /// P-value cutoff for enriched GO terms (right-tailed hypergeometric)
    #[arg(long, default_value_t = DEFAULT_GO_ENRICHMENT_P_VALUE)]
    go_enrichment_p_value: f64,

    /// Minimum positive-gene support required for a GO term
    #[arg(long, default_value_t = DEFAULT_GO_ENRICHMENT_MIN_SUPPORT)]
    go_enrichment_min_support: usize,

    /// Include direct parents of enriched GO terms in the pool
    #[arg(long, default_value_t = DEFAULT_GO_ENRICHMENT_INCLUDE_PARENTS)]
    go_enrichment_include_parents: bool,

    /// Minimum absolute fold change to keep a term (keep if fold>=T or <=1/T)
    #[arg(long, default_value_t = DEFAULT_GO_ENRICHMENT_MIN_FOLD)]
    go_enrichment_min_fold: f64,

    /// Maximum background frequency allowed for a GO term (1.0 disables)
    #[arg(long, default_value_t = DEFAULT_GO_ENRICHMENT_MAX_BG_FREQ)]
    go_enrichment_max_bg_freq: f64,

    /// Filter out GO root/namespace-level terms
    #[arg(long, default_value_t = DEFAULT_GO_ENRICHMENT_FILTER_ROOTS)]
    go_enrichment_filter_roots: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let hpo_term: TermId = args
        .hpo_term
        .parse()
        .with_context(|| format!("Invalid HPO term '{}'", args.hpo_term))?;

    let metric = parse_metric(&args.metric, args.beta)?;

    let go = go_ontology();
    let gene_set = if args.use_expanded {
        println!("Using pre-expanded GO annotations (direct + ancestors).");
        gene_set_annotations_expanded(&go)
    } else {
        println!("Using direct GO annotations (runtime traversal).");
        gene_set_annotations()
    };

    let pos_genes = count_positive_genes(&gene_set, &hpo_term);
    let total_genes = gene_set.get_gene_annotations_map().len();
    println!(
        "Dataset: {} positive / {} total genes for {}",
        pos_genes, total_genes, hpo_term
    );

    let candidate_terms = build_go_pool(
        &gene_set,
        &go,
        &hpo_term,
        &args,
        args.use_enriched_go_pool,
    );
    println!(
        "GO pool size: {} (enriched: {})",
        candidate_terms.len(),
        args.use_enriched_go_pool
    );

    let config = DecisionTreeConfig {
        metric,
        max_depth: args.max_depth,
        min_support: args.min_support,
        min_improvement: args.min_improvement,
        min_leaf_score: args.min_leaf_score,
        allow_go_negations: !args.disallow_go_negations,
        candidate_terms: Some(candidate_terms),
    };

    if args.use_expanded {
        let checker = PreexpandedSatisfactionChecker::new(&gene_set);
        run_tree(&checker, &gene_set, &hpo_term, &config)?;
    } else {
        let checker = NaiveSatisfactionChecker::new(&go, &gene_set);
        run_tree(&checker, &gene_set, &hpo_term, &config)?;
    }

    Ok(())
}

fn run_tree<C: go2hpo_genetic_algorithm::logical_formula::SatisfactionChecker>(
    checker: &C,
    gene_set: &GeneSetAnnotations,
    hpo_term: &TermId,
    config: &DecisionTreeConfig,
) -> Result<()> {
    println!(
        "Training decision tree: max_depth={}, min_support={}, min_improvement={:.2e}, min_leaf_score={:.2e}",
        config.max_depth, config.min_support, config.min_improvement, config.min_leaf_score
    );

    let result: DecisionTreeResult = decision_tree_rules(checker, hpo_term, config);

    print_overall_summary(&result, gene_set, hpo_term, config);

    let mut active_leaves: Vec<&DecisionTreeLeaf> =
        result.leaves.iter().filter(|leaf| leaf.is_active).collect();
    active_leaves.sort_by(|a, b| {
        b.summary
            .score
            .partial_cmp(&a.summary.score)
            .unwrap_or(Ordering::Equal)
    });

    if active_leaves.is_empty() {
        println!("No active leaves met the thresholds.");
        return Ok(());
    }

    println!("\nActive rules (conjunctions):");
    for (idx, leaf) in active_leaves.iter().enumerate() {
        println!(
            "  #{:<2} depth={} support={} score={:.4} precision={:.4} recall={:.4} | {}",
            idx + 1,
            leaf.depth,
            leaf.predicted_support,
            leaf.summary.score,
            leaf.summary.precision,
            leaf.summary.recall,
            leaf.conjunction
        );
    }

    Ok(())
}

fn print_overall_summary(
    result: &DecisionTreeResult,
    gene_set: &GeneSetAnnotations,
    hpo_term: &TermId,
    config: &DecisionTreeConfig,
) {
    let support_pos = count_positive_genes(gene_set, hpo_term);
    let total = gene_set.get_gene_annotations_map().len();
    let support_neg = total.saturating_sub(support_pos);

    println!(
        "\nOverall summary: score={:.4} precision={:.4} recall={:.4} | pos_support={} neg_support={}",
        result.summary.score,
        result.summary.precision,
        result.summary.recall,
        support_pos,
        support_neg
    );

    let active = result.leaves.iter().filter(|leaf| leaf.is_active).count();
    println!(
        "Leaves: {} total, {} active (min_leaf_score >= {:.2e}, min_support >= {})",
        result.leaves.len(),
        active,
        config.min_leaf_score,
        config.min_support
    );
}

fn parse_metric(raw: &str, beta: f64) -> Result<ScoreMetric> {
    let lower = raw.trim().to_lowercase();
    match lower.as_str() {
        "accuracy" => Ok(ScoreMetric::Accuracy),
        "precision" => Ok(ScoreMetric::Precision),
        "recall" => Ok(ScoreMetric::Recall),
        "f" | "fscore" | "f1" => Ok(ScoreMetric::FScore(1.0)),
        other if other.starts_with('f') => {
            let beta_str = other.trim_start_matches('f');
            let parsed_beta = if beta_str.is_empty() {
                beta
            } else {
                beta_str.parse::<f64>().with_context(|| {
                    format!("Could not parse beta from metric '{}'", raw)
                })?
            };
            Ok(ScoreMetric::FScore(parsed_beta))
        }
        _ => Err(anyhow!(
            "Unsupported metric '{}'. Use accuracy|precision|recall|f|f1|f2|f0.5|f{beta}",
            raw
        )),
    }
}

fn count_positive_genes(gene_set: &GeneSetAnnotations, phenotype: &TermId) -> usize {
    gene_set
        .get_gene_annotations_map()
        .values()
        .filter(|ann| ann.contains_phenotype(phenotype))
        .count()
}

fn build_go_pool(
    gene_set: &GeneSetAnnotations,
    go: &MinimalCsrOntology,
    hpo_term: &TermId,
    args: &Args,
    use_enriched: bool,
) -> Vec<TermId> {
    let union_terms: Vec<TermId> = gene_set
        .get_go_terms_for_hpo_phenotype(hpo_term)
        .into_iter()
        .collect();

    if !use_enriched {
        if union_terms.is_empty() {
            println!("⚠ GO pool is empty; falling back to full ontology.");
            return go.iter_term_ids().cloned().collect();
        }
        return union_terms;
    }

    let (enriched_terms, table) = compute_enriched_go_terms(
        gene_set,
        go,
        hpo_term,
        args.go_enrichment_top_k,
        args.go_enrichment_p_value,
        args.go_enrichment_min_support,
        args.go_enrichment_include_parents,
        args.go_enrichment_min_fold,
        args.go_enrichment_max_bg_freq,
        args.go_enrichment_filter_roots,
    );

    if enriched_terms.is_empty() {
        println!(
            "⚠ Enriched GO pool empty (top_k={}, min_support={}, p_value_cutoff={:.3}); falling back to union ({} terms).",
            args.go_enrichment_top_k,
            args.go_enrichment_min_support,
            args.go_enrichment_p_value,
            union_terms.len()
        );
        return if union_terms.is_empty() {
            go.iter_term_ids().cloned().collect()
        } else {
            union_terms
        };
    }

    println!(
        "Using enriched GO pool: selected {} terms (top_k={}, min_support={}, p_value_cutoff={:.3}, min_fold={:.2}, max_bg_freq={:.2}, filter_roots={}, include_parents={})",
        enriched_terms.len(),
        args.go_enrichment_top_k,
        args.go_enrichment_min_support,
        args.go_enrichment_p_value,
        args.go_enrichment_min_fold,
        args.go_enrichment_max_bg_freq,
        args.go_enrichment_filter_roots,
        args.go_enrichment_include_parents
    );

    if !table.is_empty() {
        println!("Top enriched GO terms:");
        for entry in table.iter().take(5) {
            let name = go
                .term_by_id(&entry.term)
                .map(|t| t.name().to_string())
                .unwrap_or_else(|| "(name not found)".to_string());
            println!(
                "  {} ({}) | pos={} total={} | p={:.2e} | fold={:.2}",
                entry.term, name, entry.pos_with_term, entry.total_with_term, entry.p_value, entry.fold_enrichment
            );
        }
    }

    enriched_terms
}

#[derive(Debug, Clone)]
struct GoTermEnrichment {
    term: TermId,
    p_value: f64,
    fold_enrichment: f64,
    pos_with_term: usize,
    total_with_term: usize,
}

fn compute_enriched_go_terms(
    gene_set_annotations: &GeneSetAnnotations,
    go_ontology: &MinimalCsrOntology,
    hpo_term: &TermId,
    top_k: usize,
    p_value_cutoff: f64,
    min_positive_support: usize,
    include_parents: bool,
    min_fold: f64,
    max_bg_freq: f64,
    filter_roots: bool,
) -> (Vec<TermId>, Vec<GoTermEnrichment>) {
    let gene_map = gene_set_annotations.get_gene_annotations_map();
    let total_genes = gene_map.len();
    if total_genes == 0 {
        return (Vec::new(), Vec::new());
    }

    let mut pos_total = 0usize;
    let mut total_counts: HashMap<TermId, usize> = HashMap::new();
    let mut pos_counts: HashMap<TermId, usize> = HashMap::new();

    for ann in gene_map.values() {
        let is_pos = ann.contains_phenotype(hpo_term);
        if is_pos {
            pos_total += 1;
        }

        for term in ann.get_term_annotations() {
            *total_counts.entry(term.clone()).or_default() += 1;
            if is_pos {
                *pos_counts.entry(term.clone()).or_default() += 1;
            }
        }
    }

    if pos_total == 0 {
        return (Vec::new(), Vec::new());
    }

    let neg_total = total_genes.saturating_sub(pos_total);

    let root_set: HashSet<TermId> = if filter_roots {
        [
            "GO:0008150",
            "GO:0009987",
            "GO:0003674",
            "GO:0005575",
            "GO:0110165",
        ]
        .iter()
        .filter_map(|id| id.parse::<TermId>().ok())
        .collect()
    } else {
        HashSet::new()
    };

    let mut scores: Vec<GoTermEnrichment> = total_counts
        .into_iter()
        .filter_map(|(term, total_with_term)| {
            if filter_roots && root_set.contains(&term) {
                return None;
            }

            let bg_freq = total_with_term as f64 / total_genes as f64;
            if max_bg_freq < 1.0 && bg_freq > max_bg_freq {
                return None;
            }

            let pos_with_term = *pos_counts.get(&term).unwrap_or(&0);
            if pos_with_term < min_positive_support {
                return None;
            }

            let neg_with_term = total_with_term.saturating_sub(pos_with_term);
            let pos_rate = pos_with_term as f64 / pos_total as f64;
            let neg_rate = if neg_total == 0 {
                0.0
            } else {
                neg_with_term as f64 / neg_total as f64
            };
            let fold_enrichment = if neg_rate == 0.0 {
                f64::INFINITY
            } else {
                pos_rate / neg_rate
            };

            if min_fold > 1.0 {
                let inv = 1.0 / min_fold;
                if !(fold_enrichment >= min_fold || fold_enrichment <= inv) {
                    return None;
                }
            }

            let successes = total_with_term as u64;
            let failures = (total_genes - total_with_term) as u64;
            let draws = pos_total as u64;
            let p_value = if successes == 0 || draws == 0 {
                1.0
            } else {
                statrs::distribution::Hypergeometric::new(successes, failures, draws)
                    .ok()
                    .map(|hg| 1.0 - hg.cdf((pos_with_term as u64).saturating_sub(1)))
                    .unwrap_or(1.0)
            };

            Some(GoTermEnrichment {
                term,
                p_value,
                fold_enrichment,
                pos_with_term,
                total_with_term,
            })
        })
        .collect();

    scores.sort_by(|a, b| {
        a.p_value
            .partial_cmp(&b.p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                b.fold_enrichment
                    .partial_cmp(&a.fold_enrichment)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| b.pos_with_term.cmp(&a.pos_with_term))
    });

    let mut selected: Vec<TermId> = scores
        .iter()
        .filter(|s| s.p_value <= p_value_cutoff)
        .take(top_k)
        .map(|s| s.term.clone())
        .collect();

    if selected.is_empty() {
        selected = scores
            .iter()
            .take(top_k)
            .map(|s| s.term.clone())
            .collect();
    }

    if include_parents && !selected.is_empty() {
        let mut with_parents: HashSet<TermId> = selected.iter().cloned().collect();
        for term in selected.iter() {
            for parent in go_ontology.iter_parent_ids(term) {
                with_parents.insert(parent.clone());
            }
        }
        selected = with_parents.into_iter().collect();
    }

    (selected, scores)
}

