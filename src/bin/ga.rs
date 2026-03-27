mod ga_common;

use clap::Parser;
use ga_common::{
    run_ga, GaConfig, DEFAULT_GO_ENRICHMENT_FILTER_ROOTS, DEFAULT_GO_ENRICHMENT_INCLUDE_PARENTS,
    DEFAULT_GO_ENRICHMENT_MAX_BG_FREQ, DEFAULT_GO_ENRICHMENT_MIN_FOLD,
    DEFAULT_GO_ENRICHMENT_MIN_SUPPORT, DEFAULT_GO_ENRICHMENT_P_VALUE, DEFAULT_GO_ENRICHMENT_TOP_K,
};
use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::{
    gene_set_annotations, gene_set_annotations_expanded, go_ontology, gtex_summary, phenotype2genes,
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

    /// Export last generation as bincode snapshot
    #[arg(long)]
    export_bin: Option<String>,

    /// Import snapshot to seed the initial population
    #[arg(long)]
    import_bin: Option<String>,

    /// Use pre-expanded GO annotations (direct + ancestors) for traversal-free checking
    #[arg(long, default_value_t = false)]
    use_expanded: bool,

    /// Disallow GO-term negations (NOT(GO:...)) in generated formulas
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
    let gene_set_annotations = if args.use_expanded {
        println!("Using pre-expanded GO annotations (direct + ancestors).");
        gene_set_annotations_expanded(&go_ontology)
    } else {
        println!("Using direct GO annotations (runtime traversal).");
        gene_set_annotations()
    };
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
        export_bin: args.export_bin,
        import_bin: args.import_bin,
        use_expanded: args.use_expanded,
        allow_go_negations: !args.disallow_go_negations,
        use_enriched_go_pool: args.use_enriched_go_pool,
        go_enrichment_top_k: args.go_enrichment_top_k,
        go_enrichment_p_value: args.go_enrichment_p_value,
        go_enrichment_min_support: args.go_enrichment_min_support,
        go_enrichment_include_parents: args.go_enrichment_include_parents,
        go_enrichment_min_fold: args.go_enrichment_min_fold,
        go_enrichment_max_bg_freq: args.go_enrichment_max_bg_freq,
        go_enrichment_filter_roots: args.go_enrichment_filter_roots,
    };

    // Run the GA
    let _result = run_ga(
        &config,
        &go_ontology,
        &gtex,
        &gene_set_annotations,
        &hpo2genes,
    );
}
