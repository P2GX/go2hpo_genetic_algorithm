use std::collections::HashMap;

use ontolius::TermId;

use crate::annotations::{GeneId, GeneSetAnnotations};
use crate::genetic_algorithm::ScoreMetric;

#[derive(Debug, Clone, Copy, Default)]
pub struct ConfusionCounts {
    pub tp: usize,
    pub fp: usize,
    pub tn: usize,
    pub fn_: usize,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct MetricSummary {
    pub score: f64,
    pub precision: f64,
    pub recall: f64,
    pub accuracy: f64,
    pub support_pos: usize,
    pub support_neg: usize,
}

pub fn confusion_from_predictions(
    gene_set: &GeneSetAnnotations,
    phenotype: &TermId,
    predictions: &HashMap<GeneId, bool>,
) -> ConfusionCounts {
    gene_set
        .get_gene_annotations_map()
        .iter()
        .fold(ConfusionCounts::default(), |mut acc, (gene, ann)| {
            let predicted = *predictions.get(gene).unwrap_or(&false);
            let actual = ann.contains_phenotype(phenotype);

            match (predicted, actual) {
                (true, true) => acc.tp += 1,
                (true, false) => acc.fp += 1,
                (false, true) => acc.fn_ += 1,
                (false, false) => acc.tn += 1,
            }

            acc
        })
}

fn safe_precision(counts: &ConfusionCounts) -> f64 {
    let denom = counts.tp + counts.fp;
    if denom == 0 {
        0.0
    } else {
        counts.tp as f64 / denom as f64
    }
}

fn safe_recall(counts: &ConfusionCounts) -> f64 {
    let denom = counts.tp + counts.fn_;
    if denom == 0 {
        0.0
    } else {
        counts.tp as f64 / denom as f64
    }
}

fn safe_accuracy(counts: &ConfusionCounts) -> f64 {
    let total = counts.tp + counts.fp + counts.tn + counts.fn_;
    if total == 0 {
        0.0
    } else {
        (counts.tp + counts.tn) as f64 / total as f64
    }
}

fn f_score(counts: &ConfusionCounts, beta: f64) -> f64 {
    let precision = safe_precision(counts);
    let recall = safe_recall(counts);

    if precision + recall == 0.0 {
        return 0.0;
    }

    let beta_sq = beta * beta;
    (1.0 + beta_sq) * (precision * recall) / (beta_sq * precision + recall)
}

fn score_value(counts: &ConfusionCounts, metric: &ScoreMetric) -> f64 {
    match metric {
        ScoreMetric::Accuracy => safe_accuracy(counts),
        ScoreMetric::Precision => safe_precision(counts),
        ScoreMetric::Recall => safe_recall(counts),
        ScoreMetric::FScore(beta) => f_score(counts, *beta),
    }
}

pub fn summarize(counts: &ConfusionCounts, metric: &ScoreMetric) -> MetricSummary {
    MetricSummary {
        score: score_value(counts, metric),
        precision: safe_precision(counts),
        recall: safe_recall(counts),
        accuracy: safe_accuracy(counts),
        support_pos: counts.tp + counts.fn_,
        support_neg: counts.tn + counts.fp,
    }
}

