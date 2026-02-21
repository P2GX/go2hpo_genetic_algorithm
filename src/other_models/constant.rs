use std::collections::HashMap;

use ontolius::TermId;

use crate::{
    annotations::GeneSetAnnotations,
    genetic_algorithm::ScoreMetric,
    other_models::metrics::{confusion_from_predictions, summarize, MetricSummary},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstantStrategy {
    AllPositive,
    AllNegative,
    MetricBest,
}

#[derive(Debug, Clone)]
pub struct ConstantBaselineResult {
    pub strategy_used: ConstantStrategy,
    pub summary: MetricSummary,
}

pub fn run_constant_baseline(
    gene_set: &GeneSetAnnotations,
    phenotype: &TermId,
    metric: &ScoreMetric,
    strategy: ConstantStrategy,
) -> ConstantBaselineResult {
    let positive_summary =
        summarize_predictions(gene_set, phenotype, metric, ConstantStrategy::AllPositive);
    let negative_summary =
        summarize_predictions(gene_set, phenotype, metric, ConstantStrategy::AllNegative);

    let (strategy_used, summary) = match strategy {
        ConstantStrategy::AllPositive => (ConstantStrategy::AllPositive, positive_summary),
        ConstantStrategy::AllNegative => (ConstantStrategy::AllNegative, negative_summary),
        ConstantStrategy::MetricBest => {
            if positive_summary.score > negative_summary.score
                || (positive_summary.score - negative_summary.score).abs() < f64::EPSILON
                    && positive_summary.recall >= negative_summary.recall
            {
                (ConstantStrategy::AllPositive, positive_summary)
            } else {
                (ConstantStrategy::AllNegative, negative_summary)
            }
        }
    };

    ConstantBaselineResult {
        strategy_used,
        summary,
    }
}

fn summarize_predictions(
    gene_set: &GeneSetAnnotations,
    phenotype: &TermId,
    metric: &ScoreMetric,
    strategy: ConstantStrategy,
) -> MetricSummary {
    let prediction_value = matches!(strategy, ConstantStrategy::AllPositive);
    let predictions = build_predictions(gene_set, prediction_value);
    let counts = confusion_from_predictions(gene_set, phenotype, &predictions);
    summarize(&counts, metric)
}

fn build_predictions(gene_set: &GeneSetAnnotations, value: bool) -> HashMap<String, bool> {
    gene_set
        .get_gene_annotations_map()
        .keys()
        .map(|gene| (gene.clone(), value))
        .collect()
}

