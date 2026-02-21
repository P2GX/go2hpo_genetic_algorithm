mod constant;
mod decision_tree;
mod greedy;
mod metrics;
mod single_term;

pub use constant::{run_constant_baseline, ConstantBaselineResult, ConstantStrategy};
pub use decision_tree::{
    decision_tree_rules, DecisionTreeConfig, DecisionTreeLeaf, DecisionTreeResult,
};
pub use greedy::{greedy_forward_rule, GreedyBaselineResult, GreedyConfig, GreedyStep, LiteralKind};
pub use metrics::{confusion_from_predictions, ConfusionCounts, MetricSummary, summarize};
pub use single_term::{best_single_term, SingleTermBaselineResult, SingleTermConfig};

