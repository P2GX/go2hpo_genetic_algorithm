use csv::Writer;
use go2hpo_genetic_algorithm::{
    annotations::{GeneAnnotations, GeneSetAnnotations},
    export_import::{
        export_snapshot, import_snapshot, rebuild_population_from_snapshot, GaPopulationSnapshot,
        GaRunMetadata, SerializableSolution, SNAPSHOT_SCHEMA_VERSION,
    },
    genetic_algorithm::{
        ConjunctionMutation, ConjunctionScorer, DNFScorer, DNFVecCrossover, ElitesByNumberSelector,
        FormulaEvaluator, GeneticAlgorithm, ScoreMetric, SimpleDNFVecMutation, TournamentSelection,
    },
    logical_formula::{
        DNFVec, GenePickerConjunctionGenerator, NaiveSatisfactionChecker,
        RandomConjunctionGenerator, RandomDNFVecGenerator, SatisfactionChecker, DNF,
    },
    Solution,
};
use gtex_analyzer::expression_analysis::GtexSummary;
use ontolius::ontology::csr::MinimalCsrOntology;
use ontolius::ontology::OntologyTerms;
use ontolius::term::MinimalTerm;
use ontolius::TermId;
use rand::{rngs::SmallRng, SeedableRng};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use std::sync::Arc;

/// Configuration for running the Genetic Algorithm
#[derive(Debug, Clone)]
pub struct GaConfig {
    pub hpo_term: TermId,
    pub pop_size: usize,
    pub generations: usize,
    pub mutation_rate: f64,
    pub tournament_size: usize,
    pub max_n_terms: usize,
    pub max_n_conj: usize,
    pub penalty_lambda: f64,
    pub fscore_beta: Option<f64>, // None means use estimated beta
    pub output_file: Option<String>,
    pub rng_seed: u64,
    pub export_bin: Option<String>,
    pub import_bin: Option<String>,
}

impl Default for GaConfig {
    fn default() -> Self {
        Self {
            hpo_term: "HP:0000001".parse().unwrap(),
            pop_size: 20,
            generations: 5,
            mutation_rate: 0.5,
            tournament_size: 3,
            max_n_terms: 5,
            max_n_conj: 4,
            penalty_lambda: 0.0,
            fscore_beta: None,
            output_file: None,
            rng_seed: 42,
            export_bin: None,
            import_bin: None,
        }
    }
}

/// Estimates the optimal F-score beta value based on class imbalance.
///
/// Approaches:
/// 1. **Class imbalance ratio**: beta = sqrt(negative_count / positive_count)
///    This gives higher beta for more imbalanced datasets
/// 2. **Logarithmic scaling**: beta = log10(imbalance_ratio + 1) + 1
///    More conservative, caps beta at reasonable values
/// 3. **Fixed heuristic**: Based on imbalance ratio ranges
///
/// Returns: (recommended_beta, method_used, imbalance_ratio)
pub fn estimate_fscore_beta(positive_count: usize, total_count: usize) -> (f64, &'static str, f64) {
    if positive_count == 0 {
        return (1.0, "default (no positives)", 0.0);
    }

    let negative_count = total_count - positive_count;
    let imbalance_ratio = negative_count as f64 / positive_count as f64;

    // Method: Square root of imbalance ratio (common heuristic in imbalanced learning)
    // This gives beta proportional to how imbalanced the data is
    // Formula: beta = sqrt(negative_count / positive_count)
    // Rationale: For highly imbalanced data, we need higher beta to emphasize recall
    //
    // Alternative methods (not used but documented):
    // - Logarithmic: beta = log10(imbalance_ratio + 1) + 1 (more conservative)
    // - Fixed ranges: beta = 1.0-3.0 based on imbalance thresholds
    let recommended_beta = imbalance_ratio.sqrt().clamp(1.0, 3.0);

    (recommended_beta, "sqrt(imbalance_ratio)", imbalance_ratio)
}

/// Gets the count of genes annotated with a specific HPO term
pub fn get_hpo_gene_count(
    phenotype2genes: &HashMap<TermId, HashSet<String>>,
    hpo_term: &TermId,
) -> u32 {
    phenotype2genes
        .get(hpo_term)
        .map(|genes| genes.len() as u32)
        .unwrap_or(0)
}

/// Writes GA statistics to a CSV file, preceded by metadata comments.
/// Metadata lines start with `#` and are ignored by pandas/R when reading with `comment="#"`.
pub fn write_generation_stats_to_csv(
    path: &str,
    stats_history: &[(f64, f64, f64, usize, f64, usize, f64, f64)],
    metadata: &str, // e.g. parameters description
) -> Result<(), Box<dyn std::error::Error>> {
    // Open file for writing
    let mut file = File::create(path)?;

    // Write metadata block at top (convert newlines to "# " prefixed lines)
    if !metadata.is_empty() {
        writeln!(file, "# {}", metadata.replace('\n', "\n# "))?;
    }

    // Prepare CSV writer using same file handle
    let mut wtr = Writer::from_writer(file);

    // Header
    wtr.write_record(&[
        "generation",
        "min",
        "avg",
        "max",
        "min_len",
        "avg_len",
        "max_len",
        "best_one_precision",
        "best_one_recall",
    ])?;

    // Data rows with fixed precision (5 decimals)
    for (gen, (min, avg, max, min_len, avg_len, max_len, best_one_precision, best_one_recall)) in
        stats_history.iter().enumerate()
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

/// Runs the Genetic Algorithm with the given configuration
///
/// Returns the stats history and the best solution from the last generation
pub fn run_ga(
    config: &GaConfig,
    go_ontology: &MinimalCsrOntology,
    gtex: &GtexSummary,
    gene_set_annotations: &GeneSetAnnotations,
    hpo2genes: &HashMap<TermId, HashSet<String>>,
) -> (
    Vec<(f64, f64, f64, usize, f64, usize, f64, f64)>,
    Solution<DNFVec>,
) {
    // Elite preservation: lower values (5%) promote more diversity and exploration
    // Higher values (10-20%) preserve best solutions but may cause premature convergence
    let elite_percentage = 0.05; // Recommended: 0.05 (5%) for better exploration

    let hpo_gene_count = get_hpo_gene_count(hpo2genes, &config.hpo_term);
    println!(
        "Phenotype {} has {} genes annotated",
        config.hpo_term, hpo_gene_count
    );

    // --- Get filtered GO terms for this HPO term ---
    let filtered_go_terms_set =
        gene_set_annotations.get_go_terms_for_hpo_phenotype(&config.hpo_term);
    let filtered_go_terms: Vec<TermId> = filtered_go_terms_set.into_iter().collect();
    println!(
        "Found {} GO terms annotated to genes with phenotype {}",
        filtered_go_terms.len(),
        config.hpo_term
    );

    // --- Estimate optimal F-score beta based on class imbalance ---
    let total_gene_count = gene_set_annotations.get_gene_annotations_map().len();
    let positive_count = hpo_gene_count as usize;
    let (estimated_beta, method, imbalance_ratio) =
        estimate_fscore_beta(positive_count, total_gene_count);
    println!(
        "Class imbalance: {} positives / {} total (ratio: {:.2}:1). Estimated beta: {:.2} (method: {})",
        positive_count, total_gene_count, imbalance_ratio, estimated_beta, method
    );

    // Use configured beta or estimated beta
    let fscore_beta = config.fscore_beta.unwrap_or(estimated_beta);
    println!("Using F-score beta = {:.2}", fscore_beta);

    // --- RNGs ---
    let mut rng_main = SmallRng::seed_from_u64(config.rng_seed);
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
        gene_set_annotations,
        Some(config.hpo_term.clone()),
        Some(2),
        Some(2),
    );
    let dnf_gen = RandomDNFVecGenerator::new(&mut conj_gen, 2, rng_dnf_gen);

    // --- Evaluator ---
    let checker = Arc::new(NaiveSatisfactionChecker::new(
        go_ontology,
        gene_set_annotations,
    ));
    let conj_scorer =
        ConjunctionScorer::new(Arc::clone(&checker), ScoreMetric::FScore(fscore_beta));
    let scorer = DNFScorer::new(conj_scorer, config.penalty_lambda);
    let evaluator = FormulaEvaluator::new(Box::new(scorer));

    // --- Operators & pools ---
    // Use filtered GO terms if available, otherwise use all terms from ontology
    let go_terms_pool: Vec<_> = if !filtered_go_terms.is_empty() {
        filtered_go_terms.clone()
    } else {
        go_ontology.iter_term_ids().cloned().collect()
    };

    let tissue_terms: Vec<String> = gtex
        .metadata
        .get_tissue_names()
        .into_iter()
        .cloned()
        .collect();

    // --- Optional import of an existing population snapshot ---
    let mut imported_population: Option<Vec<Solution<DNFVec>>> = None;
    if let Some(import_path) = &config.import_bin {
        println!("Importing population from {}", import_path);
        let snapshot = match import_snapshot(import_path) {
            Ok(s) => s,
            Err(e) => {
                eprintln!("⚠ Failed to read snapshot {}: {}", import_path, e);
                std::process::exit(1);
            }
        };

        let snapshot_hpo_term: TermId = match snapshot.metadata.hpo_term.parse() {
            Ok(t) => t,
            Err(e) => {
                eprintln!(
                    "⚠ Snapshot HPO term ({}) could not be parsed: {}",
                    snapshot.metadata.hpo_term, e
                );
                std::process::exit(1);
            }
        };

        if snapshot_hpo_term != config.hpo_term {
            eprintln!(
                "⚠ Snapshot HPO term ({}) does not match requested term ({})",
                snapshot_hpo_term, config.hpo_term
            );
            std::process::exit(1);
        }

        let rebuilt_population = match rebuild_population_from_snapshot(
            &snapshot,
            &evaluator,
            &config.hpo_term,
            go_ontology,
            &tissue_terms,
        ) {
            Ok(p) => p,
            Err(e) => {
                eprintln!("⚠ Failed to rebuild population from snapshot: {}", e);
                std::process::exit(1);
            }
        };

        if snapshot.metadata.pop_size != rebuilt_population.len() {
            eprintln!(
                "⚠ Snapshot pop_size ({}) differs from reconstructed size ({}); using reconstructed size.",
                snapshot.metadata.pop_size,
                rebuilt_population.len()
            );
        }

        imported_population = Some(rebuilt_population);
    }

    let selection = Box::new(TournamentSelection::new(
        config.tournament_size,
        &mut rng_selection,
    ));
    let crossover = Box::new(DNFVecCrossover::new(&mut rng_crossover));
    let mutation = Box::new(SimpleDNFVecMutation::new(
        if !filtered_go_terms.is_empty() {
            ConjunctionMutation::new_with_filtered_go_terms(
                go_ontology,
                gtex,
                config.max_n_terms,
                &mut rng_conj_mut,
                &filtered_go_terms,
            )
        } else {
            ConjunctionMutation::new(go_ontology, gtex, config.max_n_terms, &mut rng_conj_mut)
        },
        RandomConjunctionGenerator::new(1, &go_terms_pool, 1, &tissue_terms, rng_main.clone()),
        config.max_n_conj,
        &mut rng_disj_mut,
    ));

    let pop_size_for_run = imported_population
        .as_ref()
        .map(|p| p.len())
        .unwrap_or(config.pop_size);
    let numb_elites = ((pop_size_for_run as f64 * elite_percentage).ceil() as usize).max(1);

    println!(
        "\nRunning GA: HPO={}, pop_size={}, gens={}, mutation_rate={}, penalty_lambda={}, tournament_size={}, elites={}",
        config.hpo_term, pop_size_for_run, config.generations, config.mutation_rate,
        config.penalty_lambda, config.tournament_size, numb_elites
    );

    let elites = Box::new(ElitesByNumberSelector::new(numb_elites));

    // --- Build GA ---
    let mut ga = if let Some(population) = imported_population {
        GeneticAlgorithm::new_with_population(
            population,
            evaluator,
            selection,
            crossover,
            mutation,
            elites,
            Box::new(dnf_gen),
            config.mutation_rate,
            config.generations,
            rng_main,
            config.hpo_term.clone(),
        )
    } else {
        GeneticAlgorithm::new_with_size(
            pop_size_for_run,
            evaluator,
            selection,
            crossover,
            mutation,
            elites,
            Box::new(dnf_gen),
            config.mutation_rate,
            config.generations,
            rng_main,
            config.hpo_term.clone(),
        )
    };

    // --- Run GA ---
    let stats_history = ga.fit_with_stats_history();

    for (
        gen,
        (min, avg_score, max, min_len, avg_len, max_len, best_one_precision, best_one_recall),
    ) in stats_history.iter().enumerate()
    {
        println!(
            "Gen {}: min = {:.4}, avg = {:.4}, max = {:.4}, min_len = {}, avg_len = {:.4}, max_len = {}, best_one_precision = {}, best_one_recall = {}",
            gen, min, avg_score, max, min_len, avg_len, max_len, best_one_precision, best_one_recall
        );
    }

    // Compute best solution for metadata and later use
    let best_last = ga
        .get_population()
        .iter()
        .max_by(|a: &&Solution<DNFVec>, b| a.partial_cmp(b).unwrap())
        .unwrap()
        .clone();

    // Export binary snapshot of the last generation if requested
    if let Some(ref export_path) = config.export_bin {
        let snapshot = GaPopulationSnapshot {
            metadata: GaRunMetadata {
                schema_version: SNAPSHOT_SCHEMA_VERSION,
                hpo_term: config.hpo_term.to_string(),
                pop_size: ga.get_population().len(),
                generations: config.generations,
                mutation_rate: config.mutation_rate,
                tournament_size: config.tournament_size,
                max_n_terms: config.max_n_terms,
                max_n_conj: config.max_n_conj,
                penalty_lambda: config.penalty_lambda,
                fscore_beta,
                rng_seed: config.rng_seed,
                generation_index: config.generations,
                filtered_go_terms: if !filtered_go_terms.is_empty() {
                    Some(filtered_go_terms.iter().map(|t| t.to_string()).collect())
                } else {
                    None
                },
                tissue_terms_used: tissue_terms.clone(),
                best_score: best_last.get_score(),
            },
            population: ga
                .get_population()
                .iter()
                .map(SerializableSolution::from)
                .collect(),
        };

        match export_snapshot(export_path, &snapshot) {
            Ok(_) => println!("✔ Exported GA snapshot to {}", export_path),
            Err(e) => eprintln!("⚠ Failed to export snapshot to {}: {}", export_path, e),
        }
    }

    // Save results if output file is specified
    if let Some(ref file_name) = config.output_file {
        if !file_name.is_empty() {
            // Ensure the "stats" folder exists
            if let Err(e) = std::fs::create_dir_all("stats") {
                eprintln!("⚠ Failed to create 'stats' directory: {}", e);
            }

            // Build full path inside the folder
            let csv_path = format!("stats/{}.csv", file_name);

            // Get the formula string representation
            let best_formula_str = format!("{}", best_last.get_formula());

            let metadata = format!(
                "HPO term: {}\nPopulation size: {}\nGenerations: {}\nMutation rate: {}\nTournament size: {}\nMax terms per Conjunction: {}\nMax conjunctions in DNF: {}\nPenalty lambda: {}\nF-score beta: {:.2}\nBest formula: {}",
                config.hpo_term,
                pop_size_for_run,
                config.generations,
                config.mutation_rate,
                config.tournament_size,
                config.max_n_terms,
                config.max_n_conj,
                config.penalty_lambda,
                fscore_beta,
                best_formula_str,
            );

            match write_generation_stats_to_csv(&csv_path, &stats_history, &metadata) {
                Ok(_) => println!("✔ Results saved to {}", csv_path),
                Err(e) => eprintln!("⚠ Failed to write CSV: {}", e),
            }
        }
    }

    // Print detailed results
    println!("Best solution (last generation) = {}\n", best_last);

    let hpo_annot_genes: HashMap<&String, &GeneAnnotations> = checker
        .get_gene_set()
        .get_gene_annotations_map()
        .iter()
        .filter(|(_, ann)| ann.contains_phenotype(&config.hpo_term))
        .collect();

    println!(
        "Total genes annotated to {} = {}",
        config.hpo_term,
        hpo_annot_genes.len()
    );

    for conj in best_last.get_formula().get_active_conjunctions() {
        let satisfied: HashMap<String, bool> = checker.all_satisfactions(conj);
        let num_all_satisfied = satisfied.values().filter(|&&v| v).count();

        // count how many HPO-annotated genes are satisfied
        let mut num_hpo_satisfied = 0;
        for (gene_id, _) in &hpo_annot_genes {
            if checker.is_satisfied(gene_id, conj) {
                num_hpo_satisfied += 1;
            }
        }

        println!(
            "Conjunction {} satisfied by {} genes, of which {} annot. to the hpo term",
            conj, num_all_satisfied, num_hpo_satisfied
        );

        println!("Term names in conjunction:");
        for term in &conj.term_observations {
            let term_id = &term.term_id;
            if let Some(term) = go_ontology.term_by_id(term_id) {
                println!("{}: {}", term_id, term.name());
            } else {
                println!("{}: (name not found)", term_id);
            }
        }
        println!("--------------------------------");
    }

    (stats_history, best_last)
}
