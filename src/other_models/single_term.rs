use std::collections::HashSet;

use ontolius::TermId;

use crate::{
    genetic_algorithm::ScoreMetric,
    logical_formula::{Conjunction, SatisfactionChecker, TermObservation},
    other_models::metrics::{confusion_from_predictions, summarize, MetricSummary},
};

#[derive(Debug, Clone)]
pub struct SingleTermConfig {
    pub metric: ScoreMetric,
    pub min_support: usize,
    /// Minimum metric gain required to accept a candidate (guards noise).
    pub min_improvement: f64,
    /// If metric scores are within `min_improvement`, require this recall gain.
    pub min_recall_gain: f64,
}

impl Default for SingleTermConfig {
    fn default() -> Self {
        Self {
            // Emphasize recall-heavy scenarios by default.
            metric: ScoreMetric::FScore(2.5),
            min_support: 1,
            min_improvement: 1e-6,
            min_recall_gain: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct SingleTermBaselineResult {
    pub best_term: Option<TermId>,
    pub summary: MetricSummary,
}

pub fn best_single_term<C: SatisfactionChecker>(
    checker: &C,
    phenotype: &TermId,
    config: &SingleTermConfig,
) -> SingleTermBaselineResult {
    let gene_set = checker.get_gene_set();

    let candidate_terms: HashSet<TermId> = gene_set
        .get_gene_annotations_map()
        .values()
        .flat_map(|ann| ann.get_term_annotations().iter().cloned())
        .collect();

    let mut best_term: Option<TermId> = None;
    // Start from -inf so the first viable candidate is accepted.
    let mut best_summary = MetricSummary {
        score: f64::NEG_INFINITY,
        ..MetricSummary::default()
    };

    for term in candidate_terms {
        let conj = Conjunction::from(vec![TermObservation::new(term.clone(), false)], vec![]);
        let predictions = checker.all_satisfactions(&conj);
        let predicted_support = predictions.values().filter(|v| **v).count();
        if predicted_support < config.min_support {
            continue;
        }

        let counts = confusion_from_predictions(gene_set, phenotype, &predictions);
        let summary = summarize(&counts, &config.metric);

        if better(&summary, &best_summary, config) {
            best_term = Some(term);
            best_summary = summary;
        }
    }

    SingleTermBaselineResult {
        best_term,
        summary: best_summary,
    }
}

fn better(candidate: &MetricSummary, best: &MetricSummary, config: &SingleTermConfig) -> bool {
    let improves_metric = candidate.score > best.score + config.min_improvement;
    let recall_favored = (candidate.score - best.score).abs() <= config.min_improvement
        && candidate.recall > best.recall + config.min_recall_gain;

    if !improves_metric && !recall_favored {
        return false;
    }

    // If metric tie, prefer higher recall; otherwise rely on metric.
    if improves_metric {
        true
    } else {
        candidate.recall > best.recall + config.min_recall_gain
    }
}

