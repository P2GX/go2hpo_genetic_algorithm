use std::collections::HashSet;

use ontolius::TermId;

use crate::{
    genetic_algorithm::ScoreMetric,
    logical_formula::{Conjunction, Formula, SatisfactionChecker, TermObservation, TissueExpression},
    other_models::metrics::{confusion_from_predictions, summarize, MetricSummary},
};

#[derive(Debug, Clone, PartialEq)]
pub enum LiteralKind {
    Term(TermObservation),
    Tissue(TissueExpression),
}

#[derive(Debug, Clone)]
pub struct GreedyConfig {
    pub metric: ScoreMetric,
    pub max_literals: usize,
    pub min_improvement: f64,
    pub min_support: usize,
    /// Accept recall-driven gains even if metric gain is within tolerance.
    pub min_recall_gain: f64,
}

impl Default for GreedyConfig {
    fn default() -> Self {
        Self {
            // Default to recall-friendly scoring.
            metric: ScoreMetric::FScore(2.5),
            max_literals: 5,
            min_improvement: 1e-6,
            min_support: 1,
            min_recall_gain: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct GreedyStep {
    pub added: LiteralKind,
    pub summary: MetricSummary,
}

#[derive(Debug, Clone)]
pub struct GreedyBaselineResult {
    pub conjunction: Conjunction,
    pub summary: MetricSummary,
    pub steps: Vec<GreedyStep>,
}

pub fn greedy_forward_rule<C: SatisfactionChecker>(
    checker: &C,
    phenotype: &TermId,
    config: &GreedyConfig,
) -> GreedyBaselineResult {
    let gene_set = checker.get_gene_set();
    let mut current_conj = Conjunction::new();
    let current_predictions = checker.all_satisfactions(&current_conj);
    let mut current_summary = summarize(
        &confusion_from_predictions(gene_set, phenotype, &current_predictions),
        &config.metric,
    );
    let mut steps: Vec<GreedyStep> = Vec::new();

    let mut term_candidates = collect_term_candidates(gene_set);
    let mut tissue_candidates = collect_tissue_candidates(gene_set);

    while current_conj.len() < config.max_literals {
        let mut best_candidate: Option<(LiteralKind, MetricSummary, Conjunction)> = None;

        for term in term_candidates.iter() {
            if current_conj.term_observations.contains(term) {
                continue;
            }

            let mut new_conj = current_conj.clone();
            new_conj.term_observations.push(term.clone());
            let predictions = checker.all_satisfactions(&new_conj);
            let predicted_support = predictions.values().filter(|p| **p).count();
            if predicted_support < config.min_support {
                continue;
            }

            let counts = confusion_from_predictions(gene_set, phenotype, &predictions);
            let summary = summarize(&counts, &config.metric);

            if is_better(
                &summary,
                &current_summary,
                best_candidate.as_ref(),
                config.min_improvement,
                config.min_recall_gain,
            ) {
                best_candidate = Some((LiteralKind::Term(term.clone()), summary, new_conj));
            }
        }

        for tissue in tissue_candidates.iter() {
            if current_conj.tissue_expressions.contains(tissue) {
                continue;
            }

            let mut new_conj = current_conj.clone();
            new_conj.tissue_expressions.push(tissue.clone());
            let predictions = checker.all_satisfactions(&new_conj);
            let predicted_support = predictions.values().filter(|p| **p).count();
            if predicted_support < config.min_support {
                continue;
            }

            let counts = confusion_from_predictions(gene_set, phenotype, &predictions);
            let summary = summarize(&counts, &config.metric);

            if is_better(
                &summary,
                &current_summary,
                best_candidate.as_ref(),
                config.min_improvement,
                config.min_recall_gain,
            ) {
                best_candidate = Some((LiteralKind::Tissue(tissue.clone()), summary, new_conj));
            }
        }

        match best_candidate {
            Some((literal, summary, new_conj)) => {
                current_conj = new_conj;
                current_summary = summary;
                steps.push(GreedyStep {
                    added: literal,
                    summary,
                });

                // Optional: remove used candidate to avoid re-evaluating it.
                term_candidates.retain(|c| !current_conj.term_observations.contains(c));
                tissue_candidates.retain(|c| !current_conj.tissue_expressions.contains(c));
            }
            None => break,
        }
    }

    GreedyBaselineResult {
        conjunction: current_conj,
        summary: current_summary,
        steps,
    }
}

fn collect_term_candidates(gene_set: &crate::annotations::GeneSetAnnotations) -> Vec<TermObservation> {
    let mut set: HashSet<TermId> = HashSet::new();
    for ann in gene_set.get_gene_annotations_map().values() {
        set.extend(ann.get_term_annotations().iter().cloned());
    }
    set.into_iter().map(|t| TermObservation::new(t, false)).collect()
}

fn collect_tissue_candidates(
    gene_set: &crate::annotations::GeneSetAnnotations,
) -> Vec<TissueExpression> {
    let mut set: HashSet<TissueExpression> = HashSet::new();
    for ann in gene_set.get_gene_annotations_map().values() {
        set.extend(ann.get_tissue_expressions().iter().cloned());
    }
    set.into_iter().collect()
}

fn is_better(
    candidate: &MetricSummary,
    current: &MetricSummary,
    best_so_far: Option<&(LiteralKind, MetricSummary, Conjunction)>,
    min_improvement: f64,
    min_recall_gain: f64,
) -> bool {
    let improves_metric = candidate.score > current.score + min_improvement;
    let acceptable_metric = candidate.score + min_improvement >= current.score;
    let recall_gain = candidate.recall > current.recall + min_recall_gain;
    let recall_driven = acceptable_metric && recall_gain;

    if !improves_metric && !recall_driven {
        return false;
    }

    match best_so_far {
        None => true,
        Some((_, best_summary, _)) => {
            if improves_metric {
                candidate.score > best_summary.score
                    || (candidate.score - best_summary.score).abs() < f64::EPSILON
                        && candidate.recall > best_summary.recall
            } else {
                // Recall-driven branch: allow choosing the candidate with more recall
                (candidate.score + min_improvement >= best_summary.score)
                    && candidate.recall > best_summary.recall
            }
        }
    }
}

