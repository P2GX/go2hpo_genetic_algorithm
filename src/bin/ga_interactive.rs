mod ga_common;

use ga_common::{estimate_fscore_beta, run_ga, GaConfig};
use go2hpo_genetic_algorithm::annotations::GeneSetAnnotations;
use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::{go_ontology, gtex_summary, gene_set_annotations, phenotype2genes};
use gtex_analyzer::expression_analysis::GtexSummary;
use ontolius::ontology::csr::MinimalCsrOntology;
use ontolius::TermId;
use std::io::{self, Write};




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

        // Estimate beta first to show to user
        let hpo_gene_count = ga_common::get_hpo_gene_count(&hpo2genes, &hpo_term);
        let total_gene_count = gene_set_annotations.get_gene_annotations_map().len();
        let positive_count = hpo_gene_count as usize;
        let (estimated_beta, method, imbalance_ratio) = estimate_fscore_beta(positive_count, total_gene_count);
        
        println!(
            "Class imbalance: {} positives / {} total (ratio: {:.2}:1). Estimated beta: {:.2} (method: {})",
            positive_count, total_gene_count, imbalance_ratio, estimated_beta, method
        );

        // Ask user if they want to use the estimated beta or override it
        print!("Use estimated beta ({:.2})? [Y/n, or enter custom value]: ", estimated_beta);
        io::stdout().flush().unwrap();
        let mut beta_in = String::new();
        io::stdin().read_line(&mut beta_in).unwrap();
        let beta_input = beta_in.trim();
        
        let fscore_beta: Option<f64> = if beta_input.is_empty() || beta_input.eq_ignore_ascii_case("y") || beta_input.eq_ignore_ascii_case("yes") {
            None // Use estimated (will be set in run_ga)
        } else {
            match beta_input.parse() {
                Ok(beta) => Some(beta),
                Err(_) => {
                    println!("Invalid input, using estimated beta {:.2}", estimated_beta);
                    None
                }
            }
        };

        // Build configuration
        let config = GaConfig {
            hpo_term,
            pop_size,
            generations,
            mutation_rate,
            tournament_size,
            max_n_terms,
            max_n_conj,
            penalty_lambda,
            fscore_beta,
            output_file: if file_name.is_empty() { None } else { Some(file_name) },
            rng_seed: 42,
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
    println!("Exited interactive GA session.");
}


// TO RUN IT
// cargo run --bin ga_interactive