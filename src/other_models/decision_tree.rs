use std::collections::{HashMap, HashSet};

use ontolius::TermId;

use crate::{
    annotations::GeneSetAnnotations,
    genetic_algorithm::ScoreMetric,
    logical_formula::{Conjunction, SatisfactionChecker, TermObservation},
};

use super::metrics::{confusion_from_predictions, summarize, MetricSummary};

/// Configuration for the decision-tree baseline.
#[derive(Debug, Clone)]
pub struct DecisionTreeConfig {
    pub metric: ScoreMetric,
    /// Maximum depth of the tree (root at depth 0).
    pub max_depth: usize,
    /// Minimum number of genes a rule must cover to be considered.
    pub min_support: usize,
    /// Minimum improvement in score required to accept a split.
    pub min_improvement: f64,
    /// Minimum score for a leaf to be kept as an active rule.
    pub min_leaf_score: f64,
    /// Allow NOT(GO:...) literals. When false, only positive literals are used.
    pub allow_go_negations: bool,
    /// Optional pre-filtered GO term pool. If None, use union of annotated terms.
    pub candidate_terms: Option<Vec<TermId>>,
}

impl Default for DecisionTreeConfig {
    fn default() -> Self {
        Self {
            metric: ScoreMetric::FScore(1.0),
            max_depth: 4,
            min_support: 1,
            min_improvement: 1e-6,
            // Require a tiny positive score so degenerate zero-precision leaves are ignored.
            min_leaf_score: 1e-6,
            allow_go_negations: true,
            candidate_terms: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DecisionTreeLeaf {
    pub conjunction: Conjunction,
    pub summary: MetricSummary,
    pub depth: usize,
    pub predicted_support: usize,
    pub is_active: bool,
}

#[derive(Debug, Clone)]
pub struct DecisionTreeResult {
    pub leaves: Vec<DecisionTreeLeaf>,
    pub summary: MetricSummary,
}

/// Learn a small decision tree that yields a set of conjunction rules (DNF).
/// The tree only splits on GO term presence/absence to keep conjunctions representable.
pub fn decision_tree_rules<C: SatisfactionChecker>(
    checker: &C,
    phenotype: &TermId,
    config: &DecisionTreeConfig,
) -> DecisionTreeResult {
    let gene_set = checker.get_gene_set();
    let candidate_terms: Vec<TermId> = if let Some(custom) = &config.candidate_terms {
        custom.iter().cloned().collect::<HashSet<_>>().into_iter().collect()
    } else {
        gene_set
            .get_gene_annotations_map()
            .values()
            .flat_map(|ann| ann.get_term_annotations().iter())
            .cloned()
            .collect::<HashSet<_>>()
            .into_iter()
            .collect()
    };

    let mut leaves = Vec::new();
    let mut nodes_visited = 0usize;
    grow_tree(
        checker,
        gene_set,
        phenotype,
        config,
        &candidate_terms,
        Conjunction::new(),
        0,
        HashSet::new(),
        &mut leaves,
        &mut nodes_visited,
    );

    let summary = summarize_overall(checker, gene_set, phenotype, config, &leaves);

    DecisionTreeResult { leaves, summary }
}

struct SplitChoice {
    term: TermId,
    conj_true: Conjunction,
    conj_false: Option<Conjunction>,
    support_true: usize,
    support_false: usize,
    summary_true: MetricSummary,
    summary_false: Option<MetricSummary>,
}

fn grow_tree<C: SatisfactionChecker>(
    checker: &C,
    gene_set: &GeneSetAnnotations,
    phenotype: &TermId,
    config: &DecisionTreeConfig,
    candidate_terms: &[TermId],
    path: Conjunction,
    depth: usize,
    used_terms: HashSet<TermId>,
    leaves: &mut Vec<DecisionTreeLeaf>,
    nodes_visited: &mut usize,
) {
    *nodes_visited += 1;
    if *nodes_visited % 100 == 0 {
        println!("Visited {} nodes (depth {})...", nodes_visited, depth);
    }

    let predictions = checker.all_satisfactions(&path);
    let support = positive_support(&predictions);
    let counts = confusion_from_predictions(gene_set, phenotype, &predictions);
    let summary = summarize(&counts, &config.metric);

    let depth_limit_reached = depth >= config.max_depth;
    let support_too_small = support < config.min_support;
    let exhausted_terms = used_terms.len() >= candidate_terms.len();

    if depth_limit_reached || support_too_small || exhausted_terms {
        leaves.push(make_leaf(path, summary, depth, support, config));
        return;
    }

    if let Some(split) = best_split(
        checker,
        gene_set,
        phenotype,
        config,
        candidate_terms,
        &used_terms,
        &path,
        &summary,
        config.allow_go_negations,
    ) {
        let mut used_terms_next = used_terms.clone();
        used_terms_next.insert(split.term.clone());

        grow_tree(
            checker,
            gene_set,
            phenotype,
            config,
            candidate_terms,
            split.conj_true,
            depth + 1,
            used_terms_next.clone(),
            leaves,
            nodes_visited,
        );

        if let Some(conj_false) = split.conj_false {
            grow_tree(
                checker,
                gene_set,
                phenotype,
                config,
                candidate_terms,
                conj_false,
                depth + 1,
                used_terms_next,
                leaves,
                nodes_visited,
            );
        }
    } else {
        leaves.push(make_leaf(path, summary, depth, support, config));
    }
}

fn best_split<C: SatisfactionChecker>(
    checker: &C,
    gene_set: &GeneSetAnnotations,
    phenotype: &TermId,
    config: &DecisionTreeConfig,
    candidate_terms: &[TermId],
    used_terms: &HashSet<TermId>,
    path: &Conjunction,
    current_summary: &MetricSummary,
    allow_negations: bool,
) -> Option<SplitChoice> {
    let mut best_choice: Option<(SplitChoice, f64)> = None;

    for term in candidate_terms {
        if used_terms.contains(term) {
            continue;
        }

        let mut conj_true = path.clone();
        conj_true
            .term_observations
            .push(TermObservation::new(term.clone(), false));
        let preds_true = checker.all_satisfactions(&conj_true);
        let support_true = positive_support(&preds_true);
        let summary_true =
            summarize(&confusion_from_predictions(gene_set, phenotype, &preds_true), &config.metric);

        let (conj_false, support_false, summary_false) = if allow_negations {
            let mut cf = path.clone();
            cf.term_observations
                .push(TermObservation::new(term.clone(), true));
            let preds_false = checker.all_satisfactions(&cf);
            let support_false = positive_support(&preds_false);
            let summary_false = summarize(
                &confusion_from_predictions(gene_set, phenotype, &preds_false),
                &config.metric,
            );
            (Some(cf), support_false, Some(summary_false))
        } else {
            (None, 0usize, None)
        };

        let combined_summary = combined_child_summary(
            gene_set,
            phenotype,
            config,
            (&preds_true, &summary_true, support_true),
            conj_false
                .as_ref()
                .zip(summary_false.as_ref())
                .map(|(cf, sf)| {
                    let preds_false = checker.all_satisfactions(cf);
                    (preds_false, sf, support_false)
                }),
        );

        let improves = combined_summary.score > current_summary.score + config.min_improvement;

        if improves {
            match &mut best_choice {
                None => {
                    best_choice = Some((
                        SplitChoice {
                            term: term.clone(),
                            conj_true,
                            conj_false,
                            support_true,
                            support_false,
                            summary_true,
                            summary_false,
                        },
                        combined_summary.score,
                    ));
                }
                Some((_, best_score)) => {
                    if combined_summary.score > *best_score
                        || (combined_summary.score - *best_score).abs() < f64::EPSILON
                            && combined_summary.recall > current_summary.recall
                    {
                        *best_score = combined_summary.score;
                        best_choice = Some((
                            SplitChoice {
                                term: term.clone(),
                                conj_true,
                                conj_false,
                                support_true,
                                support_false,
                                summary_true,
                                summary_false,
                            },
                            combined_summary.score,
                        ));
                    }
                }
            }
        }
    }

    best_choice.map(|(choice, _)| choice)
}

fn combined_child_summary(
    gene_set: &GeneSetAnnotations,
    phenotype: &TermId,
    config: &DecisionTreeConfig,
    true_branch: (&HashMap<String, bool>, &MetricSummary, usize),
    false_branch: Option<(HashMap<String, bool>, &MetricSummary, usize)>,
) -> MetricSummary {
    let mut predictions: HashMap<String, bool> = HashMap::new();

    if is_active(true_branch.1, true_branch.2, config) {
        merge_predictions(true_branch.0, &mut predictions);
    }

    if let Some((preds_false, summary_false, support_false)) = false_branch {
        if is_active(summary_false, support_false, config) {
            merge_predictions(&preds_false, &mut predictions);
        }
    }

    let counts = confusion_from_predictions(gene_set, phenotype, &predictions);
    summarize(&counts, &config.metric)
}

fn merge_predictions(source: &HashMap<String, bool>, target: &mut HashMap<String, bool>) {
    for (gene, satisfied) in source {
        if *satisfied {
            target.insert(gene.clone(), true);
        }
    }
}

fn make_leaf(
    conjunction: Conjunction,
    summary: MetricSummary,
    depth: usize,
    support: usize,
    config: &DecisionTreeConfig,
) -> DecisionTreeLeaf {
    DecisionTreeLeaf {
        conjunction,
        summary,
        depth,
        predicted_support: support,
        is_active: is_active(&summary, support, config),
    }
}

fn is_active(summary: &MetricSummary, support: usize, config: &DecisionTreeConfig) -> bool {
    support >= config.min_support
        && support > 0
        && summary.score >= config.min_leaf_score
}

fn positive_support(predictions: &HashMap<String, bool>) -> usize {
    predictions.values().filter(|p| **p).count()
}

fn summarize_overall<C: SatisfactionChecker>(
    checker: &C,
    gene_set: &GeneSetAnnotations,
    phenotype: &TermId,
    config: &DecisionTreeConfig,
    leaves: &[DecisionTreeLeaf],
) -> MetricSummary {
    let mut predictions: HashMap<String, bool> = HashMap::new();

    for leaf in leaves.iter().filter(|leaf| leaf.is_active) {
        let satisfactions = checker.all_satisfactions(&leaf.conjunction);
        merge_predictions(&satisfactions, &mut predictions);
    }

    let counts = confusion_from_predictions(gene_set, phenotype, &predictions);
    summarize(&counts, &config.metric)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::annotations::{GeneAnnotations, GeneId};
    use crate::logical_formula::SatisfactionChecker;

    fn term(id: &str) -> TermId {
        id.parse().unwrap()
    }

    fn make_gene(
        id: &str,
        terms: &[&str],
        phenotypes: &[&str],
    ) -> (GeneId, GeneAnnotations) {
        let term_set = terms.iter().map(|t| term(t)).collect();
        let pheno_set = phenotypes.iter().map(|p| term(p)).collect();
        (
            id.to_string(),
            GeneAnnotations::new(id.to_string(), term_set, HashSet::new(), pheno_set),
        )
    }

    #[derive(Clone)]
    struct DummyChecker {
        gene_set: GeneSetAnnotations,
    }

    impl SatisfactionChecker for DummyChecker {
        fn is_satisfied(&self, symbol: &GeneId, conjunction: &Conjunction) -> bool {
            let ann = self
                .gene_set
                .get_gene_annotations(symbol)
                .expect("gene exists");

            for term_obs in &conjunction.term_observations {
                let has_term = ann.contains_term_annotation(&term_obs.term_id);
                if term_obs.is_excluded && has_term {
                    return false;
                }
                if !term_obs.is_excluded && !has_term {
                    return false;
                }
            }

            // Ignore tissues for this dummy checker
            true
        }

        fn all_satisfactions(&self, conjunction: &Conjunction) -> HashMap<String, bool> {
            self.gene_set
                .get_gene_annotations_map()
                .iter()
                .map(|(gene, _)| (gene.clone(), self.is_satisfied(gene, conjunction)))
                .collect()
        }

        fn get_gene_set(&self) -> &GeneSetAnnotations {
            &self.gene_set
        }
    }

    #[test]
    fn builds_rules_that_recover_positives() {
        let genes = HashMap::from([
            make_gene("g1", &["GO:0001"], &["HPO:0000001"]),
            make_gene("g2", &["GO:0001", "GO:0002"], &["HPO:0000001"]),
            make_gene("g3", &["GO:0002"], &[]),
            make_gene("g4", &[], &[]),
        ]);

        let gene_set = GeneSetAnnotations::new(genes);
        let checker = DummyChecker { gene_set };

        let mut config = DecisionTreeConfig::default();
        config.metric = ScoreMetric::FScore(1.0);
        config.max_depth = 3;
        config.allow_go_negations = true;

        let result = decision_tree_rules(&checker, &term("HPO:0000001"), &config);

        let active_leaves: Vec<&DecisionTreeLeaf> =
            result.leaves.iter().filter(|leaf| leaf.is_active).collect();

        assert_eq!(active_leaves.len(), 1);
        assert!(
            result.summary.score > 0.99,
            "expected near-perfect score, got {}",
            result.summary.score
        );
        assert!(
            result.summary.recall > 0.99,
            "recall should be near perfect"
        );
    }
}

