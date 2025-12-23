mod ga_common;

use std::{fs, path::Path};

use anyhow::anyhow;
use ga_common::{
    run_ga, GaConfig, GaRunResult, DEFAULT_GO_ENRICHMENT_FILTER_ROOTS,
    DEFAULT_GO_ENRICHMENT_INCLUDE_PARENTS, DEFAULT_GO_ENRICHMENT_MAX_BG_FREQ,
    DEFAULT_GO_ENRICHMENT_MIN_FOLD, DEFAULT_GO_ENRICHMENT_MIN_SUPPORT,
    DEFAULT_GO_ENRICHMENT_P_VALUE, DEFAULT_GO_ENRICHMENT_TOP_K,
};
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
    #[serde(default = "default_allow_go_negations")]
    allow_go_negations: bool,
    #[serde(default = "default_use_enriched_go_pool")]
    use_enriched_go_pool: bool,
    #[serde(default = "default_go_enrichment_top_k")]
    go_enrichment_top_k: usize,
    #[serde(default = "default_go_enrichment_p_value")]
    go_enrichment_p_value: f64,
    #[serde(default = "default_go_enrichment_min_support")]
    go_enrichment_min_support: usize,
    #[serde(default = "default_go_enrichment_include_parents")]
    go_enrichment_include_parents: bool,
    #[serde(default = "default_go_enrichment_min_fold")]
    go_enrichment_min_fold: f64,
    #[serde(default = "default_go_enrichment_max_bg_freq")]
    go_enrichment_max_bg_freq: f64,
    #[serde(default = "default_go_enrichment_filter_roots")]
    go_enrichment_filter_roots: bool,
}

#[derive(Debug, Clone)]
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
    total_hpo_genes: usize,
    satisfied_hpo_genes: usize,
    satisfied_non_hpo_genes: usize,
}

fn format_summary_table(summaries: &[RunSummary]) -> String {
    if summaries.is_empty() {
        return String::new();
    }

    let mut out = String::new();
    out.push_str("\n=== Last-generation stats summary ===\n");
    out.push_str(&format!(
        "{:<12} {:<22} {:>8} {:>8} {:>8} {:>7} {:>7} {:>7} {:>10} {:>10} {:>10} {:>10} {:>12}\n",
        "HPO term",
        "Label",
        "min",
        "avg",
        "max",
        "min_len",
        "avg_len",
        "max_len",
        "best_prec",
        "best_rec",
        "hpo_total",
        "sat_hpo",
        "sat_non_hpo"
    ));
    out.push_str(&"-".repeat(
        12 + 1 + 22 + 1 + 8 * 3 + 1 + 7 * 3 + 1 + 10 * 2 + 1 + 10 + 1 + 10 + 1 + 12,
    ));
    out.push('\n');

    for s in summaries {
        out.push_str(&format!(
            "{:<12} {:<22} {:>8.4} {:>8.4} {:>8.4} {:>7} {:>7.2} {:>7} {:>10.4} {:>10.4} {:>10} {:>10} {:>12}\n",
            s.hpo_term,
            s.label.chars().take(20).collect::<String>(),
            s.min,
            s.avg,
            s.max,
            s.min_len,
            s.avg_len,
            s.max_len,
            s.best_one_precision,
            s.best_one_recall,
            s.total_hpo_genes,
            s.satisfied_hpo_genes,
            s.satisfied_non_hpo_genes
        ));
    }

    out
}

fn default_rng_seed() -> u64 {
    42
}

fn default_allow_go_negations() -> bool {
    true
}

fn default_use_enriched_go_pool() -> bool {
    false
}

fn default_go_enrichment_top_k() -> usize {
    DEFAULT_GO_ENRICHMENT_TOP_K
}

fn default_go_enrichment_p_value() -> f64 {
    DEFAULT_GO_ENRICHMENT_P_VALUE
}

fn default_go_enrichment_min_support() -> usize {
    DEFAULT_GO_ENRICHMENT_MIN_SUPPORT
}

fn default_go_enrichment_include_parents() -> bool {
    DEFAULT_GO_ENRICHMENT_INCLUDE_PARENTS
}

fn default_go_enrichment_min_fold() -> f64 {
    DEFAULT_GO_ENRICHMENT_MIN_FOLD
}

fn default_go_enrichment_max_bg_freq() -> f64 {
    DEFAULT_GO_ENRICHMENT_MAX_BG_FREQ
}

fn default_go_enrichment_filter_roots() -> bool {
    DEFAULT_GO_ENRICHMENT_FILTER_ROOTS
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
            allow_go_negations: self.allow_go_negations,
            use_enriched_go_pool: self.use_enriched_go_pool,
            go_enrichment_top_k: self.go_enrichment_top_k,
            go_enrichment_p_value: self.go_enrichment_p_value,
            go_enrichment_min_support: self.go_enrichment_min_support,
            go_enrichment_include_parents: self.go_enrichment_include_parents,
            go_enrichment_min_fold: self.go_enrichment_min_fold,
            go_enrichment_max_bg_freq: self.go_enrichment_max_bg_freq,
            go_enrichment_filter_roots: self.go_enrichment_filter_roots,
        })
    }
}

fn resolve_summary_path(raw: &str) -> String {
    let path = Path::new(raw);
    if path.is_relative() {
        format!("stats/runs/summaries/{}", raw)
    } else {
        raw.to_string()
    }
}

fn main() -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();

    let mut config_path: Option<String> = None;
    let mut summary_path: Option<String> = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-s" | "--summary-file" => {
                if i + 1 >= args.len() {
                    anyhow::bail!("--summary-file requires a path argument");
                }
                summary_path = Some(resolve_summary_path(&args[i + 1]));
                i += 2;
            }
            other if config_path.is_none() => {
                config_path = Some(other.to_string());
                i += 1;
            }
            other => {
                anyhow::bail!(
                    "Unrecognized argument '{}'. Usage: ga_batch [config_path] [--summary-file <path>]",
                    other
                );
            }
        }
    }

    let config_path = config_path.unwrap_or_else(|| "src/bin/run_configs/hpo_runs.ljson".to_string());

    println!("Reading run list from {}", config_path);
    let config_path_for_err = config_path.clone();
    let config_text = fs::read_to_string(&config_path)
        .map_err(|e| anyhow!("Failed to read {}: {}", config_path_for_err, e))?;

    let config_path_for_err = config_path.clone();
    let runs: Vec<BatchRun> = serde_json::from_str(&config_text)
        .map_err(|e| anyhow!("Failed to parse run config {}: {}", config_path_for_err, e))?;

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
        let run_result: GaRunResult = run_ga(
            &config,
            &go_ontology,
            &gtex,
            &gene_set_annotations,
            &hpo2genes,
        );

        if let Some((min, avg, max, min_len, avg_len, max_len, best_one_precision, best_one_recall)) =
            run_result.stats_history.last().copied()
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
                total_hpo_genes: run_result.total_hpo_genes,
                satisfied_hpo_genes: run_result.satisfied_hpo_genes,
                satisfied_non_hpo_genes: run_result.satisfied_non_hpo_genes,
            });
        } else {
            println!("Run {} had no stats history; skipping summary row.", label);
        }
    }

    if !summaries.is_empty() {
        let table = format_summary_table(&summaries);
        print!("{}", table);

        if let Some(path) = summary_path {
            if let Some(parent) = Path::new(&path).parent() {
                if let Err(e) = fs::create_dir_all(parent) {
                    eprintln!("⚠ Failed to create summary directory {:?}: {}", parent, e);
                }
            }
            match fs::write(&path, &table) {
                Ok(_) => println!("Summary saved to {}", path),
                Err(e) => eprintln!("⚠ Failed to write summary to {}: {}", path, e),
            }
        }
    }

    println!("\nAll runs completed.");
    Ok(())
}

