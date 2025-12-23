mod ga_common;

use ga_common::{
    estimate_fscore_beta, run_ga, GaConfig, DEFAULT_GO_ENRICHMENT_FILTER_ROOTS,
    DEFAULT_GO_ENRICHMENT_INCLUDE_PARENTS, DEFAULT_GO_ENRICHMENT_MAX_BG_FREQ,
    DEFAULT_GO_ENRICHMENT_MIN_FOLD, DEFAULT_GO_ENRICHMENT_MIN_SUPPORT,
    DEFAULT_GO_ENRICHMENT_P_VALUE, DEFAULT_GO_ENRICHMENT_TOP_K,
};
use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::{
    gene_set_annotations, gene_set_annotations_expanded, go_ontology, gtex_summary, phenotype2genes,
};
use gtex_analyzer::expression_analysis::GtexSummary;
use ontolius::ontology::csr::MinimalCsrOntology;
use ontolius::TermId;
use std::io::{self, Write};

fn main() {
    let go_ontology: MinimalCsrOntology = go_ontology();
    let gtex: GtexSummary = gtex_summary().expect("GTEx summary should load correctly");
    println!("Use pre-expanded GO annotations? [y/N]: ");
    io::stdout().flush().unwrap();
    let mut expanded_in = String::new();
    io::stdin().read_line(&mut expanded_in).unwrap();
    let use_expanded = matches!(
        expanded_in.trim().to_lowercase().as_str(),
        "y" | "yes" | "true" | "1"
    );
    println!("Disallow GO-term negations (NOT(GO:...))? [y/N]: ");
    io::stdout().flush().unwrap();
    let mut neg_in = String::new();
    io::stdin().read_line(&mut neg_in).unwrap();
    let allow_go_negations = !matches!(
        neg_in.trim().to_lowercase().as_str(),
        "y" | "yes" | "true" | "1"
    );

    println!(
        "Use enriched GO pool (top {} by p<={:.3}, min_support={})? [y/N]: ",
        DEFAULT_GO_ENRICHMENT_TOP_K, DEFAULT_GO_ENRICHMENT_P_VALUE, DEFAULT_GO_ENRICHMENT_MIN_SUPPORT
    );
    io::stdout().flush().unwrap();
    let mut enriched_in = String::new();
    io::stdin().read_line(&mut enriched_in).unwrap();
    let use_enriched_go_pool = matches!(
        enriched_in.trim().to_lowercase().as_str(),
        "y" | "yes" | "true" | "1"
    );

    let (
        go_enrichment_top_k,
        go_enrichment_p_value,
        go_enrichment_min_support,
        go_enrichment_include_parents,
        go_enrichment_min_fold,
        go_enrichment_max_bg_freq,
        go_enrichment_filter_roots,
    ) = if use_enriched_go_pool {
        print!(
            "Enter max enriched GO terms [default={}]: ",
            DEFAULT_GO_ENRICHMENT_TOP_K
        );
        io::stdout().flush().unwrap();
        let mut topk_in = String::new();
        io::stdin().read_line(&mut topk_in).unwrap();
        let go_enrichment_top_k: usize = topk_in
            .trim()
            .parse()
            .unwrap_or(DEFAULT_GO_ENRICHMENT_TOP_K);

        print!(
            "Enter p-value cutoff (right-tailed hypergeometric) [default={:.3}]: ",
            DEFAULT_GO_ENRICHMENT_P_VALUE
        );
        io::stdout().flush().unwrap();
        let mut pval_in = String::new();
        io::stdin().read_line(&mut pval_in).unwrap();
        let go_enrichment_p_value: f64 = pval_in
            .trim()
            .parse()
            .unwrap_or(DEFAULT_GO_ENRICHMENT_P_VALUE);

        print!(
            "Enter minimum positive-gene support [default={}]: ",
            DEFAULT_GO_ENRICHMENT_MIN_SUPPORT
        );
        io::stdout().flush().unwrap();
        let mut minsupp_in = String::new();
        io::stdin().read_line(&mut minsupp_in).unwrap();
        let go_enrichment_min_support: usize = minsupp_in
            .trim()
            .parse()
            .unwrap_or(DEFAULT_GO_ENRICHMENT_MIN_SUPPORT);

        print!(
            "Enter minimum fold (keep if fold>=X or <=1/X) [default={:.2}]: ",
            DEFAULT_GO_ENRICHMENT_MIN_FOLD
        );
        io::stdout().flush().unwrap();
        let mut fold_in = String::new();
        io::stdin().read_line(&mut fold_in).unwrap();
        let go_enrichment_min_fold: f64 = fold_in
            .trim()
            .parse()
            .unwrap_or(DEFAULT_GO_ENRICHMENT_MIN_FOLD);

        print!(
            "Include parents of enriched GO terms? [Y/n, default={}]: ",
            if DEFAULT_GO_ENRICHMENT_INCLUDE_PARENTS {
                "Y"
            } else {
                "n"
            }
        );
        io::stdout().flush().unwrap();
        let mut parents_in = String::new();
        io::stdin().read_line(&mut parents_in).unwrap();
        let go_enrichment_include_parents = match parents_in.trim().to_lowercase().as_str() {
            "" => DEFAULT_GO_ENRICHMENT_INCLUDE_PARENTS,
            "y" | "yes" | "true" | "1" => true,
            "n" | "no" | "false" | "0" => false,
            _ => DEFAULT_GO_ENRICHMENT_INCLUDE_PARENTS,
        };

        print!(
            "Enter maximum background frequency to keep a term (1.0 disables) [default={:.2}]: ",
            DEFAULT_GO_ENRICHMENT_MAX_BG_FREQ
        );
        io::stdout().flush().unwrap();
        let mut bg_in = String::new();
        io::stdin().read_line(&mut bg_in).unwrap();
        let go_enrichment_max_bg_freq: f64 = bg_in
            .trim()
            .parse()
            .unwrap_or(DEFAULT_GO_ENRICHMENT_MAX_BG_FREQ);

        print!(
            "Filter GO root/namespace terms? [Y/n, default={}]: ",
            if DEFAULT_GO_ENRICHMENT_FILTER_ROOTS {
                "Y"
            } else {
                "n"
            }
        );
        io::stdout().flush().unwrap();
        let mut roots_in = String::new();
        io::stdin().read_line(&mut roots_in).unwrap();
        let go_enrichment_filter_roots = match roots_in.trim().to_lowercase().as_str() {
            "" => DEFAULT_GO_ENRICHMENT_FILTER_ROOTS,
            "y" | "yes" | "true" | "1" => true,
            "n" | "no" | "false" | "0" => false,
            _ => DEFAULT_GO_ENRICHMENT_FILTER_ROOTS,
        };

        (
            go_enrichment_top_k,
            go_enrichment_p_value,
            go_enrichment_min_support,
            go_enrichment_include_parents,
            go_enrichment_min_fold,
            go_enrichment_max_bg_freq,
            go_enrichment_filter_roots,
        )
    } else {
        (
            DEFAULT_GO_ENRICHMENT_TOP_K,
            DEFAULT_GO_ENRICHMENT_P_VALUE,
            DEFAULT_GO_ENRICHMENT_MIN_SUPPORT,
            DEFAULT_GO_ENRICHMENT_INCLUDE_PARENTS,
            DEFAULT_GO_ENRICHMENT_MIN_FOLD,
            DEFAULT_GO_ENRICHMENT_MAX_BG_FREQ,
            DEFAULT_GO_ENRICHMENT_FILTER_ROOTS,
        )
    };

    let gene_set_annotations = if use_expanded {
        println!("Using pre-expanded GO annotations (direct + ancestors).");
        gene_set_annotations_expanded(&go_ontology)
    } else {
        println!("Using direct GO annotations (runtime traversal).");
        gene_set_annotations()
    };
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
        let (estimated_beta, method, imbalance_ratio) =
            estimate_fscore_beta(positive_count, total_gene_count);

        println!(
            "Class imbalance: {} positives / {} total (ratio: {:.2}:1). Estimated beta: {:.2} (method: {})",
            positive_count, total_gene_count, imbalance_ratio, estimated_beta, method
        );

        // Ask user if they want to use the estimated beta or override it
        print!(
            "Use estimated beta ({:.2})? [Y/n, or enter custom value]: ",
            estimated_beta
        );
        io::stdout().flush().unwrap();
        let mut beta_in = String::new();
        io::stdin().read_line(&mut beta_in).unwrap();
        let beta_input = beta_in.trim();

        let fscore_beta: Option<f64> = if beta_input.is_empty()
            || beta_input.eq_ignore_ascii_case("y")
            || beta_input.eq_ignore_ascii_case("yes")
        {
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
            output_file: if file_name.is_empty() {
                None
            } else {
                Some(file_name)
            },
            rng_seed: 42,
            export_bin: None,
            import_bin: None,
            use_expanded,
            allow_go_negations,
            use_enriched_go_pool,
            go_enrichment_top_k,
            go_enrichment_p_value,
            go_enrichment_min_support,
            go_enrichment_include_parents,
            go_enrichment_min_fold,
            go_enrichment_max_bg_freq,
            go_enrichment_filter_roots,
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
    println!("Exited interactive GA session.");
}

// TO RUN IT
// cargo run --bin ga_interactive
