mod ga_common;

use clap::Parser;
use ga_common::{run_ga, GaConfig};
use go2hpo_genetic_algorithm::annotations::GeneSetAnnotations;
use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::{
    go_ontology, gtex_summary, gene_set_annotations, phenotype2genes,
};
use gtex_analyzer::expression_analysis::GtexSummary;
use ontolius::ontology::csr::MinimalCsrOntology;
use ontolius::TermId;

/// Genetic Algorithm for finding DNF logic formulas that predict HPO term annotations
#[derive(Parser, Debug)]
#[command(name = "ga")]
#[command(about = "Run genetic algorithm to find DNF formulas for HPO term prediction", long_about = None)]
struct Args {
    /// HPO term to analyze (e.g., HP:0001083)
    #[arg(long)]
    hpo_term: String,

    /// Population size
    #[arg(short = 'p', long, default_value_t = 20)]
    pop_size: usize,

    /// Number of generations
    #[arg(short = 'g', long, default_value_t = 5)]
    generations: usize,

    /// Mutation rate (0.0-1.0)
    #[arg(short = 'm', long, default_value_t = 0.5)]
    mutation_rate: f64,

    /// Tournament size for selection
    #[arg(short = 't', long, default_value_t = 3)]
    tournament_size: usize,

    /// Maximum number of terms per conjunction
    #[arg(long, default_value_t = 5)]
    max_n_terms: usize,

    /// Maximum number of conjunctions per DNF
    #[arg(long, default_value_t = 4)]
    max_n_conj: usize,

    /// Penalty lambda on DNF length
    #[arg(short = 'l', long, default_value_t = 0.0)]
    penalty_lambda: f64,

    /// F-score beta value (if not specified, will be estimated from class imbalance)
    #[arg(short = 'b', long)]
    fscore_beta: Option<f64>,

    /// Output file name (without extension, saved to stats/ directory)
    #[arg(short = 'o', long)]
    output_file: Option<String>,

    /// Random seed for reproducibility
    #[arg(long, default_value_t = 42)]
    rng_seed: u64,
}

fn main() {
    let args = Args::parse();

    // Parse HPO term
    let hpo_term: TermId = match args.hpo_term.parse() {
        Ok(term) => term,
        Err(e) => {
            eprintln!("Error: Invalid HPO term '{}': {}", args.hpo_term, e);
            std::process::exit(1);
        }
    };

    // Load data
    println!("Loading data...");
    let go_ontology: MinimalCsrOntology = go_ontology();
    let gtex: GtexSummary = gtex_summary().expect("GTEx summary should load correctly");
    let gene_set_annotations: GeneSetAnnotations = gene_set_annotations();
    let hpo2genes = phenotype2genes();
    println!("Data loaded successfully.\n");

    // Build configuration
    let config = GaConfig {
        hpo_term,
        pop_size: args.pop_size,
        generations: args.generations,
        mutation_rate: args.mutation_rate,
        tournament_size: args.tournament_size,
        max_n_terms: args.max_n_terms,
        max_n_conj: args.max_n_conj,
        penalty_lambda: args.penalty_lambda,
        fscore_beta: args.fscore_beta,
        output_file: args.output_file,
        rng_seed: args.rng_seed,
    };

    // Run the GA
    let (_stats_history, _best_solution) = run_ga(
        &config,
        &go_ontology,
        &gtex,
        &gene_set_annotations,
        &hpo2genes,
    );
}

