use go2hpo_genetic_algorithm::{
    annotations::GeneSetAnnotations,
    genetic_algorithm::{
        ConjunctionMutation, ConjunctionScorer, DNFScorer, DNFVecCrossover, ElitesByNumberSelector, FormulaEvaluator, GeneticAlgorithm, Mutation, ScoreMetric, Selection, SimpleDNFVecMutation, TournamentSelection
    },
    logical_formula::{
        Conjunction, GenePickerConjunctionGenerator, NaiveSatisfactionChecker, RandomConjunctionGenerator, RandomDNFVecGenerator
    }, utils::fixtures::gene_set_annotations::phenotype2genes,
};
use gtex_analyzer::expression_analysis::GtexSummary;

use ontolius::ontology::csr::MinimalCsrOntology;
use ontolius::ontology::OntologyTerms;
use ontolius::TermId;

use rand::{rngs::SmallRng, SeedableRng};
use std::{collections::{HashMap, HashSet}, io::{self, Write}};

use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::{
    go_ontology, gtex_summary, gene_set_annotations,
};



fn get_hpo_gene_count(
    phenotype2genes: &HashMap<TermId, HashSet<String>>,
    hpo_term: &TermId,
) -> u32 {
    phenotype2genes
        .get(hpo_term)              // Option<&HashSet<String>>
        .map(|genes| genes.len() as u32)  // convert length to u64
        .unwrap_or(0)               // return 0 if not found
}


fn main() {
    // --- Load data once ---
    let go_ontology: MinimalCsrOntology = go_ontology();
    let gtex: GtexSummary = gtex_summary().expect("GTEx summary should load correctly");
    let gene_set_annotations: GeneSetAnnotations = gene_set_annotations();
    let hpo2genes = phenotype2genes();

    println!("Data loaded. You can now run the GA on multiple HPO terms.");
    println!("Type 'quit' at any time to stop.\n");

    loop {
        // --- Ask interactively ---
        print!("Enter HPO term (e.g., HP:0001083): ");
        io::stdout().flush().unwrap();
        let mut hpo_input = String::new();
        io::stdin().read_line(&mut hpo_input).unwrap();
        let hpo_input = hpo_input.trim();
        if hpo_input.eq_ignore_ascii_case("quit") {
            break;
        }
        let hpo_term: TermId = match hpo_input.parse() {
            Ok(t) => t,
            Err(_) => {
                println!("Invalid HPO term, try again.\n");
                continue;
            }
        };

        // Population size
        print!("Enter population size [default=20]: ");
        io::stdout().flush().unwrap();
        let mut pop_in = String::new();
        io::stdin().read_line(&mut pop_in).unwrap();
        let pop_size: usize = pop_in.trim().parse().unwrap_or(20);

        // Generations
        print!("Enter number of generations [default=5]: ");
        io::stdout().flush().unwrap();
        let mut gens_in = String::new();
        io::stdin().read_line(&mut gens_in).unwrap();
        let generations: usize = gens_in.trim().parse().unwrap_or(5);

        // Mutation rate
        print!("Enter mutation rate (0.0–1.0) [default=0.2]: ");
        io::stdout().flush().unwrap();
        let mut mut_in = String::new();
        io::stdin().read_line(&mut mut_in).unwrap();
        let mutation_rate: f64 = mut_in.trim().parse().unwrap_or(0.2);

        // Penalty lambda
        print!("Enter penalty lambda [default=0.1]: ");
        io::stdout().flush().unwrap();
        let mut pen_in = String::new();
        io::stdin().read_line(&mut pen_in).unwrap();
        let penalty_lambda: f64 = pen_in.trim().parse().unwrap_or(0.1);

    
        println!(
            "\nRunning GA: HPO={}, pop_size={}, gens={}, mutation_rate={}",
            hpo_term, pop_size, generations, mutation_rate
        );

        let hpo_gene_count = get_hpo_gene_count(&hpo2genes, &hpo_term);
        println!("Phenotype {} has {} genes annotated", hpo_term, hpo_gene_count);

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
        let scorer = DNFScorer::new(conj_scorer, penalty_lambda);
        let evaluator = FormulaEvaluator::new(Box::new(scorer));

        // --- Operators ---
        let go_terms: Vec<_> = go_ontology.iter_term_ids().cloned().collect();
        let tissue_terms: Vec<String> = gtex
            .metadata
            .get_tissue_names()
            .into_iter()
            .cloned()
            .collect();

        let selection = Box::new(TournamentSelection::new(2, &mut rng_selection));
        let crossover = Box::new(DNFVecCrossover::new(&mut rng_crossover));
        let mutation = Box::new(SimpleDNFVecMutation::new(
            ConjunctionMutation::new(&go_ontology, &gtex, &mut rng_conj_mut),
            RandomConjunctionGenerator::new(1, &go_terms, 1, &tissue_terms, rng_main.clone()),
            5,
            &mut rng_disj_mut,
        ));
        let elites = Box::new(ElitesByNumberSelector::new(1));

        // --- Build GA ---
        let mut ga = GeneticAlgorithm::new_with_size(
            pop_size,
            evaluator,
            selection,
            crossover,
            mutation,
            elites,
            Box::new(dnf_gen),
            mutation_rate,
            generations,
            rng_main,
            hpo_term.clone(),
        );

        // --- Run GA ---
        let stats_history = ga.fit_with_stats_history();

        for (gen, (min, avg_score, max, min_len, avg_len, max_len)) in stats_history.iter().enumerate() {
            println!(
                "Gen {}: min = {:.4}, avg = {:.4}, max = {:.4}, min_len = {}, avg_len = {:.4}, max_len = {}",
                gen, min, avg_score, max, min_len, avg_len, max_len
            );
        }

        let best_last = ga
            .get_population()
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        println!("Best solution (last generation) = {}\n", best_last);
    }

    println!("Exited interactive GA session.");
}


// TO RUN IT
// cargo run --bin ga_interactive