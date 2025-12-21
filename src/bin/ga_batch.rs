mod ga_common;

use std::fs;

use anyhow::anyhow;
use ga_common::{run_ga, GaConfig};
use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::{
    gene_set_annotations_expanded, go_ontology, gtex_summary, phenotype2genes,
};
use ontolius::TermId;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct BatchRun {
    /// Optional human-friendly label for logs.
    #[serde(default)]
    label: Option<String>,
    hpo_term: String,
    pop_size: usize,
    generations: usize,
    mutation_rate: f64,
    tournament_size: usize,
    max_n_terms: usize,
    max_n_conj: usize,
    penalty_lambda: f64,
    #[serde(default)]
    fscore_beta: Option<f64>,
    #[serde(default)]
    output_file: Option<String>,
    #[serde(default = "default_rng_seed")]
    rng_seed: u64,
}

fn default_rng_seed() -> u64 {
    42
}

impl BatchRun {
    fn to_config(&self) -> anyhow::Result<GaConfig> {
        let hpo_term: TermId = self
            .hpo_term
            .parse()
            .map_err(|e| anyhow!("Invalid HPO term '{}': {}", self.hpo_term, e))?;

        Ok(GaConfig {
            hpo_term,
            pop_size: self.pop_size,
            generations: self.generations,
            mutation_rate: self.mutation_rate,
            tournament_size: self.tournament_size,
            max_n_terms: self.max_n_terms,
            max_n_conj: self.max_n_conj,
            penalty_lambda: self.penalty_lambda,
            fscore_beta: self.fscore_beta,
            output_file: self.output_file.clone(),
            rng_seed: self.rng_seed,
            export_bin: None,
            import_bin: None,
            use_expanded: true, // Always run with expanded annotations for batch mode.
        })
    }
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let config_path = args
        .get(1)
        .map(|s| s.as_str())
        .unwrap_or("src/bin/run_configs/hpo_runs.ljson");

    println!("Reading run list from {}", config_path);
    let config_text =
        fs::read_to_string(config_path).map_err(|e| anyhow!("Failed to read {}: {}", config_path, e))?;

    let runs: Vec<BatchRun> = serde_json::from_str(&config_text)
        .map_err(|e| anyhow!("Failed to parse run config {}: {}", config_path, e))?;

    if runs.is_empty() {
        println!("No runs found in config. Exiting.");
        return Ok(());
    }

    println!("Loading shared data once (expanded GO annotations)...");
    let go_ontology = go_ontology();
    let gtex = gtex_summary().expect("GTEx summary should load correctly");
    let gene_set_annotations = gene_set_annotations_expanded(&go_ontology);
    let hpo2genes = phenotype2genes();

    println!("Loaded {} runs.", runs.len());

    #[derive(Debug)]
    struct RunSummary {
        label: String,
        hpo_term: String,
        min: f64,
        avg: f64,
        max: f64,
        min_len: usize,
        avg_len: f64,
        max_len: usize,
        best_one_precision: f64,
        best_one_recall: f64,
    }

    let mut summaries: Vec<RunSummary> = Vec::with_capacity(runs.len());

    for (idx, run) in runs.iter().enumerate() {
        let label = run.label.as_deref().unwrap_or("unnamed");
        println!(
            "\n=== Run {}/{}: {} ({}) ===",
            idx + 1,
            runs.len(),
            run.hpo_term,
            label
        );

        let config = run.to_config()?;

        // Execute the GA; results are printed and optionally written per run.
        let (stats_history, _best_solution) = run_ga(
            &config,
            &go_ontology,
            &gtex,
            &gene_set_annotations,
            &hpo2genes,
        );

        if let Some((min, avg, max, min_len, avg_len, max_len, best_one_precision, best_one_recall)) =
            stats_history.last().copied()
        {
            summaries.push(RunSummary {
                label: label.to_string(),
                hpo_term: run.hpo_term.clone(),
                min,
                avg,
                max,
                min_len,
                avg_len,
                max_len,
                best_one_precision,
                best_one_recall,
            });
        } else {
            println!("Run {} had no stats history; skipping summary row.", label);
        }
    }

    if !summaries.is_empty() {
        println!("\n=== Last-generation stats summary ===");
        println!(
            "{:<12} {:<22} {:>8} {:>8} {:>8} {:>7} {:>7} {:>7} {:>10} {:>10}",
            "HPO term",
            "Label",
            "min",
            "avg",
            "max",
            "min_len",
            "avg_len",
            "max_len",
            "best_prec",
            "best_rec"
        );
        println!("{}", "-".repeat(12 + 1 + 22 + 1 + 8 * 3 + 1 + 7 * 3 + 1 + 10 * 2 + 8));
        for s in summaries {
            println!(
                "{:<12} {:<22} {:>8.4} {:>8.4} {:>8.4} {:>7} {:>7.2} {:>7} {:>10.4} {:>10.4}",
                s.hpo_term,
                s.label.chars().take(20).collect::<String>(),
                s.min,
                s.avg,
                s.max,
                s.min_len,
                s.avg_len,
                s.max_len,
                s.best_one_precision,
                s.best_one_recall
            );
        }
    }

    println!("\nAll runs completed.");
    Ok(())
}

