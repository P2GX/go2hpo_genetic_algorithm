use clap::Parser;
use go2hpo_genetic_algorithm::{
    annotations::GeneSetAnnotations,
    genetic_algorithm::{
        ConjunctionMutation, ConjunctionScorer, DNFScorer, DNFVecCrossover,
        ElitesByNumberSelector, FormulaEvaluator, GeneticAlgorithm, Mutation,
        ScoreMetric, Selection, TournamentSelection, SimpleDNFVecMutation,
    },
    logical_formula::{
        Conjunction, GenePickerConjunctionGenerator, RandomConjunctionGenerator,
        RandomDNFVecGenerator, NaiveSatisfactionChecker,
    },
    utils::fixtures::gene_set_annotations::{
        go_ontology, gtex_summary, gene_set_annotations,
    },
    Solution,
};
use gtex_analyzer::expression_analysis::GtexSummary;
use ontolius::ontology::OntologyTerms;
use ontolius::ontology::csr::MinimalCsrOntology;
use ontolius::TermId;
use rand::{rngs::SmallRng, SeedableRng};

/// Simple CLI for running the GO2HPO genetic algorithm

#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Cli {
    /// Target HPO term (e.g. HP:0001083)
    #[arg(short = 't', long)]
    hpo: String,

    /// Population size
    #[arg(short = 'p', long, default_value_t = 50)]
    pop_size: usize,

    /// Number of generations
    #[arg(short = 'g', long, default_value_t = 10)]
    gens: usize,

    /// Mutation rate
    #[arg(short = 'm', long, default_value_t = 0.2)]
    mutation_rate: f64,
}


// Execution example:
// cargo run -- -t HP:0001083 -p 5 -g 10
// cargo run -- -t HP:0001250 -p 100 -g 20 -m 0.3

fn main() {
    // --- Parse CLI ---
    let cli = Cli::parse();
    let hpo_term: TermId = cli.hpo.parse().expect("Invalid HPO term format");

    println!(
        "Running GA for phenotype: {}, pop_size = {}, gens = {}, mutation_rate = {}",
        hpo_term, cli.pop_size, cli.gens, cli.mutation_rate
    );

    // --- Load data (same as before) ---
    let go_ontology: MinimalCsrOntology = go_ontology();
    let gtex = gtex_summary().expect("Tissue expression features should load correctly");
    let gene_set_annotations = gene_set_annotations();

    // --- RNGs ---
    let mut rng_main = SmallRng::seed_from_u64(42);
    let mut rng_conj_gen = rng_main.clone();
    let mut rng_dnf_gen = rng_main.clone();
    let mut rng_selection = rng_main.clone();
    let mut rng_crossover = rng_main.clone();
    let mut rng_conj_mut = rng_main.clone();
    let mut rng_disj_mut = rng_main.clone();

    // --- Generator chain ---
    let mut conj_gen = GenePickerConjunctionGenerator::new(
        &mut rng_conj_gen,
        0.5,
        0.5,
        &gene_set_annotations,
        Some(hpo_term.clone()),
        Some(2),
        Some(2),
    );
    let dnf_gen = RandomDNFVecGenerator::new(&mut conj_gen, 2, rng_dnf_gen);

    // --- Evaluator ---
    let checker = NaiveSatisfactionChecker::new(&go_ontology, &gene_set_annotations);
    let conj_scorer = ConjunctionScorer::new(checker, ScoreMetric::FScore(1.0));
    let scorer = DNFScorer::new(conj_scorer, 0.03);
    let evaluator = FormulaEvaluator::new(Box::new(scorer));

    // --- Operators ---
    let go_terms: Vec<_> = go_ontology.iter_term_ids().take(5).cloned().collect();
    let tissue_terms: Vec<String> = gtex
        .metadata
        .get_tissue_names()
        .into_iter()
        .cloned()
        .take(5)
        .collect();

    let selection = Box::new(TournamentSelection::new(2, &mut rng_selection));
    let crossover = Box::new(DNFVecCrossover::new(&mut rng_crossover));
    let mutation = Box::new(go2hpo_genetic_algorithm::genetic_algorithm::SimpleDNFVecMutation::new(
        ConjunctionMutation::new(&go_ontology, &gtex, 4, &mut rng_conj_mut),
        RandomConjunctionGenerator::new(1, &go_terms, 1, &tissue_terms, rng_main.clone()),
        6,
        &mut rng_disj_mut,
    ));
    let elites = Box::new(ElitesByNumberSelector::new(1));

    // --- Build GA ---
    let mut ga = GeneticAlgorithm::new_with_size(
        cli.pop_size,
        evaluator,
        selection,
        crossover,
        mutation,
        elites,
        Box::new(dnf_gen),
        cli.mutation_rate,
        cli.gens,
        rng_main,
        hpo_term.clone(),
    );

    // --- Run GA ---
    let stats_history = ga.fit_with_stats_history();

    for (gen, (min, avg_score, max, min_len, avg_len, max_len, max_precision, max_recall)) in stats_history.iter().enumerate() {
            println!(
                "Gen {}: min = {:.4}, avg = {:.4}, max = {:.4}, min_len = {}, avg_len = {:.4}, max_len = {}, max_precision = {}, max_recall = {}",
                gen, min, avg_score, max, min_len, avg_len, max_len, max_precision, max_recall
            );
    }

    let best_last = ga
        .get_population()
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    println!("Best solution (last generation) = {}", best_last);
}

