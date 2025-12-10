use go2hpo_genetic_algorithm::{
    annotations::{GeneAnnotations, GeneSetAnnotations},
    genetic_algorithm::{
        ConjunctionCrossover, ConjunctionMutation, ConjunctionScorer, ElitesByNumberSelector,
        FormulaEvaluator, GeneticAlgorithm, Mutation, ScoreMetric, Selection, TournamentSelection,
    },
    logical_formula::{
        Conjunction, GenePickerConjunctionGenerator, NaiveSatisfactionChecker, SatisfactionChecker,
    },
    utils::fixtures::gene_set_annotations::{
        gene_set_annotations, go_ontology, gtex_summary, phenotype2genes,
    },
    Solution,
};
use gtex_analyzer::expression_analysis::GtexSummary;

use ontolius::ontology::csr::MinimalCsrOntology;
use ontolius::TermId;

use rand::{rngs::SmallRng, SeedableRng};
use std::{
    collections::{HashMap, HashSet},
    io::{self, Write},
    sync::Arc,
};

fn get_hpo_gene_count(
    phenotype2genes: &HashMap<TermId, HashSet<String>>,
    hpo_term: &TermId,
) -> u32 {
    phenotype2genes
        .get(hpo_term)
        .map(|genes| genes.len() as u32)
        .unwrap_or(0)
}

fn main() {
    // --- Load data once ---
    let go_ontology: MinimalCsrOntology = go_ontology();
    let gtex: GtexSummary = gtex_summary().expect("GTEx summary should load correctly");
    let gene_set_annotations: GeneSetAnnotations = gene_set_annotations();
    let hpo2genes = phenotype2genes();

    println!("Data loaded. You can now run the Conjunction-only GA on multiple HPO terms.");
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

        println!(
            "\nRunning Conjunction GA: HPO={}, pop_size={}, gens={}, mutation_rate={}",
            hpo_term, pop_size, generations, mutation_rate
        );

        let hpo_gene_count = get_hpo_gene_count(&hpo2genes, &hpo_term);
        println!(
            "Phenotype {} has {} genes annotated",
            hpo_term, hpo_gene_count
        );

        // --- RNGs ---
        let mut rng_main = SmallRng::seed_from_u64(42);
        let mut rng_conj_gen = rng_main.clone();
        let mut rng_selection = rng_main.clone();
        let mut rng_crossover = rng_main.clone();
        let mut rng_conj_mut = rng_main.clone();

        // --- Generator ---
        let mut conj_gen = GenePickerConjunctionGenerator::new(
            &mut rng_conj_gen,
            0.5,
            0.5,
            &gene_set_annotations,
            Some(hpo_term.clone()),
            Some(2),
            Some(2),
        );

        // --- Evaluator ---
        let checker = Arc::new(NaiveSatisfactionChecker::new(
            &go_ontology,
            &gene_set_annotations,
        ));
        let conj_scorer = ConjunctionScorer::new(Arc::clone(&checker), ScoreMetric::FScore(2.0)); // F2 = recall-oriented
        let evaluator = FormulaEvaluator::new(Box::new(conj_scorer));

        // --- Operators ---
        let selection = Box::new(TournamentSelection::new(2, &mut rng_selection));
        let crossover = Box::new(ConjunctionCrossover::new(&mut rng_crossover));
        let mutation = Box::new(ConjunctionMutation::new(
            &go_ontology,
            &gtex,
            6,
            &mut rng_conj_mut,
        ));

        // Elites
        let numb_elites = (pop_size as f64 * 0.1).ceil() as usize;
        let elites = Box::new(ElitesByNumberSelector::new(numb_elites));

        // --- Build GA (Conjunction-only) ---
        let mut ga = GeneticAlgorithm::new_with_size(
            pop_size,
            evaluator,
            selection,
            crossover,
            mutation,
            elites,
            Box::new(conj_gen),
            mutation_rate,
            generations,
            rng_main,
            hpo_term.clone(),
        );

        // --- Run GA ---
        let stats_history = ga.fit_with_stats_history();

        for (
            gen,
            (
                min,
                avg_score,
                max,
                _min_len,
                _avg_len,
                _max_len,
                best_one_precision,
                best_one_recall,
            ),
        ) in stats_history.iter().enumerate()
        {
            println!(
                "Gen {}: min = {:.4}, avg = {:.4}, max = {:.4}, best_one_precision = {:.4}, best_one_recall = {:.4}",
                gen, min, avg_score, max, best_one_precision, best_one_recall
            );
        }

        // --- Best solution ---
        let best_last = ga
            .get_population()
            .iter()
            .max_by(|a: &&Solution<Conjunction>, b| a.partial_cmp(b).unwrap())
            .unwrap();

        println!("Best solution (last generation) = {}\n", best_last);

        // --- Check how many HPO-annotated genes are satisfied ---
        let hpo_annot_genes: HashMap<&String, &GeneAnnotations> = checker
            .get_gene_set()
            .get_gene_annotations_map()
            .iter()
            .filter(|(_, ann)| ann.contains_phenotype(&hpo_term))
            .collect();

        println!(
            "Total genes annotated to {} = {}",
            hpo_term,
            hpo_annot_genes.len()
        );

        let mut num_hpo_satisfied = 0;
        let mut num_all_satisfied = 0;

        for (gene_id, _) in checker.get_gene_set().get_gene_annotations_map().iter() {
            if checker.is_satisfied(gene_id, best_last.get_formula()) {
                num_all_satisfied += 1;
                if hpo_annot_genes.contains_key(gene_id) {
                    num_hpo_satisfied += 1;
                }
            }
        }

        println!(
            "Best conjunction satisfied {} genes in total, of which {} are annotated to {}",
            num_all_satisfied, num_hpo_satisfied, hpo_term
        );
    }

    println!("Exited interactive Conjunction GA session.");
}
