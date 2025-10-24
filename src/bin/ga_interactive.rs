use go2hpo_genetic_algorithm::{
    annotations::{GeneAnnotations, GeneSetAnnotations},
    genetic_algorithm::{
        ConjunctionMutation, ConjunctionScorer, DNFScorer, DNFVecCrossover, ElitesByNumberSelector, FitnessScorer, FormulaEvaluator, GeneticAlgorithm, Mutation, ScoreMetric, Selection, SimpleDNFVecMutation, TournamentSelection
    },
    logical_formula::{
        Conjunction, DNFVec, GenePickerConjunctionGenerator, NaiveSatisfactionChecker, RandomConjunctionGenerator, RandomDNFVecGenerator, SatisfactionChecker, DNF
    }, utils::fixtures::gene_set_annotations::phenotype2genes, Solution,
};
use gtex_analyzer::expression_analysis::GtexSummary;

use ontolius::ontology::csr::MinimalCsrOntology;
use ontolius::ontology::OntologyTerms;
use ontolius::TermId;

use rand::{rngs::SmallRng, SeedableRng};
use std::{collections::{HashMap, HashSet}, fs, io::{self, Write}};

use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::{
    go_ontology, gtex_summary, gene_set_annotations,
};

use csv::Writer;
use std::fs::File;


use std::sync::Arc;

fn get_hpo_gene_count(
    phenotype2genes: &HashMap<TermId, HashSet<String>>,
    hpo_term: &TermId,
) -> u32 {
    phenotype2genes
        .get(hpo_term)              
        .map(|genes| genes.len() as u32)  
        .unwrap_or(0)               
}


fn write_generation_stats_to_csv(
    path: &str,
    stats_history: &[(f64, f64, f64, usize, f64, usize, f64, f64)],
) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(path)?;
    let mut wtr = Writer::from_writer(file);

    // Write header
    wtr.write_record(&[
        "generation", "min", "avg", "max", "min_len", "avg_len", "max_len",
        "best_one_precision", "best_one_recall"
    ])?;

    for (gen, (min, avg, max, min_len, avg_len, max_len, best_one_precision, best_one_recall))
        in stats_history.iter().enumerate()
    {
        wtr.write_record(&[
        gen.to_string(),
        format!("{:.5}", min),
        format!("{:.5}", avg),
        format!("{:.5}", max),
        min_len.to_string(),
        format!("{:.2}", avg_len),
        max_len.to_string(),
        format!("{:.5}", best_one_precision),
        format!("{:.5}", best_one_recall),
    ])?;
    }

    wtr.flush()?;
    Ok(())
}


fn main() {
    let go_ontology: MinimalCsrOntology = go_ontology();
    let gtex: GtexSummary = gtex_summary().expect("GTEx summary should load correctly");
    let gene_set_annotations: GeneSetAnnotations = gene_set_annotations();
    let hpo2genes = phenotype2genes();

    println!("Data loaded. You can now run the GA on multiple HPO terms.");
    println!("Type 'quit' at any time to stop.\n");

    loop {
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
        print!("Enter mutation rate (0.0–1.0) [default=0.5]: ");
        io::stdout().flush().unwrap();
        let mut mut_in = String::new();
        io::stdin().read_line(&mut mut_in).unwrap();
        let mutation_rate: f64 = mut_in.trim().parse().unwrap_or(0.5);

        // Tournament size
        print!("Enter tournament size [default=3]: ");
        io::stdout().flush().unwrap();
        let mut tourn_in = String::new();
        io::stdin().read_line(&mut tourn_in).unwrap();
        let tournament_size: usize = tourn_in.trim().parse().unwrap_or(3);

        // Max number of terms in each conjunction
        print!("Enter max number of terms per conjunction [default=5]: ");
        io::stdout().flush().unwrap();
        let mut terms_in = String::new();
        io::stdin().read_line(&mut terms_in).unwrap();
        let max_n_terms: usize = terms_in.trim().parse().unwrap_or(5);

        // Max number of conjunctions per DNF
        print!("Enter max number of conjunctions per DNF [default=4]: ");
        io::stdout().flush().unwrap();
        let mut conj_in = String::new();
        io::stdin().read_line(&mut conj_in).unwrap();
        let max_n_conj: usize = conj_in.trim().parse().unwrap_or(4);

        // Penalty lambda
        print!("Enter penalty lambda on DNF length [default=0.0]: ");
        io::stdout().flush().unwrap();
        let mut pen_in = String::new();
        io::stdin().read_line(&mut pen_in).unwrap();
        let penalty_lambda: f64 = pen_in.trim().parse().unwrap_or(0.0);

        print!("Enter output file name (without extension) [press Enter to skip saving]: ");
        io::stdout().flush().unwrap();
        let mut file_in = String::new();
        io::stdin().read_line(&mut file_in).unwrap();
        let file_name = file_in.trim().to_string();

    
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
        let checker = Arc::new(NaiveSatisfactionChecker::new(&go_ontology, &gene_set_annotations));
        let conj_scorer = ConjunctionScorer::new(Arc::clone(&checker), ScoreMetric::FScore(2.0)); //2.0 to prioritize recall over precision
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

        let selection = Box::new(TournamentSelection::new(
            tournament_size, 
            &mut rng_selection));
        let crossover = Box::new(DNFVecCrossover::new(&mut rng_crossover));
        let mutation = Box::new(SimpleDNFVecMutation::new(
            ConjunctionMutation::new(&go_ontology, &gtex, max_n_terms, &mut rng_conj_mut),
            RandomConjunctionGenerator::new(1, &go_terms, 1, &tissue_terms, rng_main.clone()),
            max_n_conj,
            &mut rng_disj_mut,
        ));

        let numb_elites = (pop_size as f64 * 0.1).ceil() as usize;
        let elites = Box::new(ElitesByNumberSelector::new(numb_elites));

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

        for (gen, (min, avg_score, max, min_len, avg_len, max_len, best_one_precision, best_one_recall)) in stats_history.iter().enumerate() {
            println!(
                "Gen {}: min = {:.4}, avg = {:.4}, max = {:.4}, min_len = {}, avg_len = {:.4}, max_len = {}, best_one_precision = {}, best_one_recall = {}",
                gen, min, avg_score, max, min_len, avg_len, max_len, best_one_precision, best_one_recall
            );
        }

        if !file_name.is_empty() {
            // Ensure the "stats" folder exists
            if let Err(e) = fs::create_dir_all("stats") {
                eprintln!("⚠ Failed to create 'stats' directory: {}", e);
            }

            // Build full path inside the folder
            let csv_path = format!("stats/{}.csv", file_name);

            match write_generation_stats_to_csv(&csv_path, &stats_history) {
                Ok(_) => println!("✔ Results saved to {}", csv_path),
                Err(e) => eprintln!("⚠ Failed to write CSV: {}", e),
            }
        } else {
            println!("⚠ No file name provided — results not saved.");
        }


        let best_last = ga
            .get_population()
            .iter()
            .max_by(|a: &&Solution<DNFVec>, b| a.partial_cmp(b).unwrap())
            .unwrap();

        //Check its precision and recall as well
        let best_formula = best_last.get_formula();
        let phenotype = hpo_term.clone();

        println!("Best solution (last generation) = {}\n", best_last);
        
        let hpo_annot_genes: HashMap<&String, &GeneAnnotations> = checker.get_gene_set().get_gene_annotations_map()
                .iter()
                .filter(|(_, ann)| ann.contains_phenotype(&hpo_term))
                .collect();

        println!(
            "Total genes annotated to {} = {}",
            hpo_term,
            hpo_annot_genes.len()
        );

        for conj in best_last.get_formula().get_active_conjunctions(){
            let satisfied: HashMap<String, bool> = checker.all_satisfactions(conj);
            let num_all_satisfied = satisfied.values().filter(|&&v| v).count();

            // count how many HPO-annotated genes are satisfied
            let mut num_hpo_satisfied = 0;
            for (gene_id, ann) in &hpo_annot_genes {
                if checker.is_satisfied(gene_id, conj) {
                    num_hpo_satisfied += 1;
                }
            }

            println!("Conjunction {} satisfied by {} genes, of which {} annot. to the hpo term", conj, num_all_satisfied, num_hpo_satisfied);
        }


    }
    println!("Exited interactive GA session.");
}


// TO RUN IT
// cargo run --bin ga_interactive