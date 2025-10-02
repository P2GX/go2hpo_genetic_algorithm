use go2hpo_genetic_algorithm::annotations::GeneSetAnnotations;
use ontolius::ontology::OntologyTerms;
use ontolius::TermId;
use rand::{rng, Rng};
use rand::rngs::ThreadRng;

use rstest::rstest;

use go2hpo_genetic_algorithm::genetic_algorithm::{ConjunctionMutation, ConjunctionScorer, DNFScorer, FormulaEvaluator, GeneticAlgorithm, Mutation, ScoreMetric, Selection, SimpleDNFBitmaskMutation, SimpleDNFVecMutation};
use go2hpo_genetic_algorithm::genetic_algorithm::TournamentSelection;
use go2hpo_genetic_algorithm::genetic_algorithm::Crossover;
use go2hpo_genetic_algorithm::genetic_algorithm::{DNFVecCrossover, ConjunctionCrossover};
use go2hpo_genetic_algorithm::genetic_algorithm::ElitesSelector;
use go2hpo_genetic_algorithm::genetic_algorithm::ElitesByNumberSelector;

use go2hpo_genetic_algorithm::Solution;

use go2hpo_genetic_algorithm::logical_formula::{
    Conjunction, DNFBitmask, DNFVec, DgeState, NaiveSatisfactionChecker, RandomConjunctionGenerator, TermObservation, TissueExpression, DNF
};

use rand::{rngs::SmallRng, SeedableRng};
use gtex_analyzer::expression_analysis::GtexSummary;
use ontolius::ontology::csr::MinimalCsrOntology;


use go2hpo_genetic_algorithm::logical_formula::FormulaGenerator;
use go2hpo_genetic_algorithm::logical_formula::{GenePickerConjunctionGenerator, RandomDNFVecGenerator};

use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::{
    go_ontology,
    gtex_summary,
    gene_set_annotations,
};


#[rstest]
fn test_genetic_algorithm_sanity(
    go_ontology: MinimalCsrOntology,
    gene_set_annotations: GeneSetAnnotations,
    gtex_summary: std::io::Result<GtexSummary>,
) {
    let gtex = gtex_summary.expect("Fixture GTEx summary must be ok");

    //RNGs
    // let mut rng = SmallRng::seed_from_u64(7);
    let mut rng_main = SmallRng::seed_from_u64(7);

    // Separate RNG clones for each operator
    let mut rng_conj_gen   = rng_main.clone();
    let mut rng_dnf_gen    = rng_main.clone();
    let mut rng_selection  = rng_main.clone();
    let mut rng_crossover  = rng_main.clone();
    let mut rng_conj_mutation   = rng_main.clone();
    let mut rng_disj_mutation   = rng_main.clone();


    // Simple generator chain
    let mut conj_gen = GenePickerConjunctionGenerator::new(
        &mut rng_conj_gen,
        0.5,
        0.5,
        &gene_set_annotations,
        Some("HP:0000001".parse().unwrap()),
        Some(2),
        Some(2),
    );
    let  dnf_gen = RandomDNFVecGenerator::new(&mut conj_gen, 2, rng_dnf_gen);

    // Evaluator
    let checker = NaiveSatisfactionChecker::new(&go_ontology, &gene_set_annotations);
    let conj_scorer = ConjunctionScorer::new(checker, ScoreMetric::Accuracy);
    let scorer = DNFScorer::new(conj_scorer, 0.0);
    let evaluator = FormulaEvaluator::new(Box::new(scorer));

    // Operators (keep simple)
    let go_terms: Vec<_> = go_ontology.iter_term_ids().take(5).cloned().collect();
    let tissue_terms: Vec<String> = gtex.metadata.get_tissue_names().into_iter().cloned().take(5).collect();
    let selection = Box::new(TournamentSelection::new(2, &mut rng_selection));
    let crossover = Box::new(DNFVecCrossover::new(&mut rng_crossover));
    let mutation = Box::new(SimpleDNFVecMutation::new(
        ConjunctionMutation::new(&go_ontology, &gtex, &mut rng_conj_mutation),
        RandomConjunctionGenerator::new(
            1,
            &go_terms,
            1,
            &tissue_terms,
            rng_main.clone(),
        ),
        4,
        &mut rng_disj_mutation,
    ));
    let elites = Box::new(ElitesByNumberSelector::new(1));

    // Build GA
    let pop_size = 4;
    let mut ga = GeneticAlgorithm::new_with_size(
        pop_size,
        evaluator,
        selection,
        crossover,
        mutation,
        elites,
        Box::new(dnf_gen),
        0.2,    // mutation rate
        2, 
        rng_main,
        "HP:0000001".parse().unwrap(),  
    );

    let best = ga.fit();
    println!("Sanity GA best solution: {}", best);

    // Assertions: population size stays the same, score is finite
    assert_eq!(ga.get_population().len(), pop_size);
    assert!(best.get_score().is_finite());
}


#[rstest]
fn test_genetic_algorithm_history(
    go_ontology: MinimalCsrOntology,
    gene_set_annotations: GeneSetAnnotations,
    gtex_summary: std::io::Result<GtexSummary>,
) {
    let gtex = gtex_summary.expect("Fixture GTEx summary must be ok");

    //RNGs
    // let mut rng = SmallRng::seed_from_u64(7);
    let mut rng_main = SmallRng::seed_from_u64(7);

    // Separate RNG clones for each operator
    let mut rng_conj_gen   = rng_main.clone();
    let mut rng_dnf_gen    = rng_main.clone();
    let mut rng_selection  = rng_main.clone();
    let mut rng_crossover  = rng_main.clone();
    let mut rng_conj_mutation   = rng_main.clone();
    let mut rng_disj_mutation   = rng_main.clone();

    // Simple generator chain
    let mut conj_gen = GenePickerConjunctionGenerator::new(
        &mut rng_conj_gen,
        0.5,
        0.5,
        &gene_set_annotations,
        Some("HP:0000001".parse().unwrap()),
        Some(2),
        Some(2),
    );
    let  dnf_gen = RandomDNFVecGenerator::new(&mut conj_gen, 2, rng_dnf_gen);

    // Evaluator
    let checker = NaiveSatisfactionChecker::new(&go_ontology, &gene_set_annotations);
    let conj_scorer = ConjunctionScorer::new(checker, ScoreMetric::Accuracy);
    let scorer = DNFScorer::new(conj_scorer, 0.0);
    let evaluator = FormulaEvaluator::new(Box::new(scorer));

    // Operators (keep simple)
    let go_terms: Vec<_> = go_ontology.iter_term_ids().take(5).cloned().collect();
    let tissue_terms: Vec<String> = gtex.metadata.get_tissue_names().into_iter().cloned().take(5).collect();
    let selection = Box::new(TournamentSelection::new(2, &mut rng_selection));
    let crossover = Box::new(DNFVecCrossover::new(&mut rng_crossover));
    let mutation = Box::new(SimpleDNFVecMutation::new(
        ConjunctionMutation::new(&go_ontology, &gtex, &mut rng_conj_mutation),
        RandomConjunctionGenerator::new(
            1,
            &go_terms,
            1,
            &tissue_terms,
            rng_main.clone(),
        ),
        4,
        &mut rng_disj_mutation,
    ));
    let elites = Box::new(ElitesByNumberSelector::new(1));

    // Build GA
    let pop_size = 4;
    let mut ga = GeneticAlgorithm::new_with_size(
        pop_size,
        evaluator,
        selection,
        crossover,
        mutation,
        elites,
        Box::new(dnf_gen),
        0.2,    // mutation rate
        10, 
        rng_main,
        "HP:0000001".parse().unwrap(),  // dummy phenotype
    );

    let history = ga.fit_with_history();

    // Print best solution per generation
    for (gen, pop) in history.iter().enumerate() {
        let best = pop.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        println!("Gen {}: best solution = {}", gen, best);
    }

    // Assertions: population size stays the same, score is finite
    assert_eq!(ga.get_population().len(), pop_size);
}




#[rstest]
fn test_genetic_algorithm_stats_history(
    go_ontology: MinimalCsrOntology,
    gene_set_annotations: GeneSetAnnotations,
    gtex_summary: std::io::Result<GtexSummary>,
) {
    let gtex = gtex_summary.expect("Fixture GTEx summary must be ok");

    // HPO term
    // HP:0000001 V
    // HP:0001083
    let hpo_term: TermId = "HP:0001083".parse().unwrap();

    // RNGs
    let mut rng_main = SmallRng::seed_from_u64(7);
    let mut rng_conj_gen   = rng_main.clone();
    let mut rng_dnf_gen    = rng_main.clone();
    let mut rng_selection  = rng_main.clone();
    let mut rng_crossover  = rng_main.clone();
    let mut rng_conj_mutation   = rng_main.clone();
    let mut rng_disj_mutation   = rng_main.clone();

    // Simple generator chain
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

    // Evaluator
    let checker = NaiveSatisfactionChecker::new(&go_ontology, &gene_set_annotations);
    let conj_scorer = ConjunctionScorer::new(checker, ScoreMetric::FScore(1.0));
    let scorer = DNFScorer::new(conj_scorer, 0.0);
    let evaluator = FormulaEvaluator::new(Box::new(scorer));

    // Operators
    let go_terms: Vec<_> = go_ontology.iter_term_ids().take(5).cloned().collect();
    let tissue_terms: Vec<String> = gtex.metadata.get_tissue_names().into_iter().cloned().take(5).collect();
    let selection = Box::new(TournamentSelection::new(2, &mut rng_selection));
    let crossover = Box::new(DNFVecCrossover::new(&mut rng_crossover));
    let mutation = Box::new(SimpleDNFVecMutation::new(
        ConjunctionMutation::new(&go_ontology, &gtex, &mut rng_conj_mutation),
        RandomConjunctionGenerator::new(
            1,
            &go_terms,
            1,
            &tissue_terms,
            rng_main.clone(),
        ),
        4,
        &mut rng_disj_mutation,
    ));
    let elites = Box::new(ElitesByNumberSelector::new(1));

    // Build GA
    let pop_size = 1;
    let mut ga = GeneticAlgorithm::new_with_size(
        pop_size,
        evaluator,
        selection,
        crossover,
        mutation,
        elites,
        Box::new(dnf_gen),
        0.2,    // mutation rate
        1,     // generations
        rng_main,
        hpo_term, // pass it here
    );

    let stats_history = ga.fit_with_stats_history();

    // Print stats for each generation
    for (gen, (min, avg_score, max, min_len, avg_len, max_len, max_precision, max_recall)) in stats_history.iter().enumerate() {
            println!(
                "Gen {}: min = {:.4}, avg = {:.4}, max = {:.4}, min_len = {}, avg_len = {:.4}, max_len = {}, max_precision = {}, max_recall = {}",
                gen, min, avg_score, max, min_len, avg_len, max_len, max_precision, max_recall
            );
        }

    // Print the best solution from the last generation
    let best_last = ga
        .get_population()
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    println!("Best solution (last generation) = {}", best_last);

    // Assertions: population size stays the same at the end
    assert_eq!(ga.get_population().len(), pop_size);
}















// #[rstest]
// fn test_genetic_algorithm_end_to_end(
//     go_ontology: MinimalCsrOntology,
//     gene_set_annotations: GeneSetAnnotations,
//     gtex_summary_sample: std::io::Result<GtexSummary>,
// ) {
//     let gtex = gtex_summary_sample.expect("Fixture GTEx summary must be ok");
//     let mut rng = SmallRng::seed_from_u64(2025);
//     let mut rng2 = SmallRng::seed_from_u64(2025);
//     let mut rng3 = SmallRng::seed_from_u64(2025);

//     // 1) Build generator chain: GenePicker -> DNFVec
//     let target_hpo: TermId = "HP:0001250".parse().unwrap(); // Example phenotype
//     let mut conj_gen = GenePickerConjunctionGenerator::new(
//         &mut rng,
//         0.7,              // prob_terms
//         0.7,              // prob_tissues
//         &gene_set_annotations,
//         Some(target_hpo.clone()), // restrict to phenotype
//         Some(3),          // max terms
//         Some(3),          // max tissues
//     );

//     let terms_mutation_pool: Vec<TermId> = go_ontology.iter_term_ids().take(10).cloned().collect();
//     let tissues_mutation_pool: Vec<String> = gtex.metadata.get_tissue_names().into_iter().take(10).cloned().collect();

//     let mut dnf_gen = RandomDNFVecGenerator::new(&mut conj_gen, 3, rng.clone());

//     // 2) Build evaluator (ConjunctionScorer -> DNFScorer -> FormulaEvaluator)
//     let checker = NaiveSatisfactionChecker::new(go_ontology, &gene_set_annotations);
//     let conj_scorer = ConjunctionScorer::new(checker, ScoreMetric::Accuracy);
//     let scorer = DNFScorer::new(conj_scorer);
//     let evaluator = FormulaEvaluator::new(Box::new(scorer));

//     // 3) Operators
//     let selection = Box::new(TournamentSelection::new(2, &mut rng2));
//     let crossover = Box::new(DNFVecCrossover::new(&mut rng3));
//     let mutation = Box::new(SimpleDNFVecMutation::new(
//         ConjunctionMutation::new(&go_ontology, &gtex, &mut rng),
//         RandomConjunctionGenerator::new(
//             2,
//             &terms_mutation_pool,
//             2,
//             &tissues_mutation_pool,
//             rng.clone(),
//         ),
//         &mut rng,
//     ));
//     let elites = Box::new(ElitesByNumberSelector::new(1));

//     // 4) Build GA with new_with_size
//     let mut ga = GeneticAlgorithm::new_with_size(
//         6,              // population size
//         evaluator,
//         selection,
//         crossover,
//         mutation,
//         elites,
//         Box::new(dnf_gen),
//         0.3,            // mutation rate
//         5,              // generations
//         rng,
//         target_hpo,     // phenotype
//     );

//     // 5) Run
//     let best = ga.fit();
//     println!("Best solution: {}", best);

//     // 6) Assert the result is valid
//     assert!(best.get_score().is_finite());
// }















// ==========================
// 5) (Optional) GA stub
// ==========================

// This is a minimal smoke test template. Uncomment and fill once you confirm
// the exact `GeneticAlgorithm::new(...)` signature and the simplest Mutation
// type you want to plug (e.g., a SimpleDNFVecMutation + a generator).


// use go2hpo_genetic_algorithm::genetic_algorithm::GeneticAlgorithm;
// use go2hpo_genetic_algorithm::SimpleDNFVecMutation;
// use go2hpo_genetic_algorithm::logical_formula::formula_generator::{
//     RedundantRandomConjunctionGenerator, RandomDNFVecGenerator, FormulaGenerator,
// };

// #[rstest]
// #[ignore] // remove this once the constructor signature is confirmed
// fn ga_end_to_end_smoke() {
//     // Fixtures (adjust to your actual fixture helpers if needed)
//     let go = fixtures::gene_set_annotations::load_and_get_go();
//     let gtex = fixtures::gene_set_annotations::load_gtex_summary();
//     let gene_set = fixtures::gene_set_annotations::GENE_SET; // or similar

//     let mut rng: ThreadRng = rng();

//     // 1) simplest generator chain (pick tiny sizes)
//     let mut conj_gen = RedundantRandomConjunctionGenerator::new(&gene_set, &go, &gtex, 1, 0, &mut rng);
//     let dnf_gen  = RandomDNFVecGenerator::new(&mut conj_gen, 1, &mut rng);

//     // // 2) evaluator (FormulaEvaluator over Conjunction/ DNFs) and scorer:
//     // let checker = ... // your SatisfactionChecker from fixtures
//     let conj_scorer = ConjunctionScorer::new(checker, &gene_set, ScoreMetric::Accuracy);
//     let scorer = DNFScorer { conjunction_scorer: conj_scorer };
//     let evaluator = FormulaEvaluator { scorer: Box::new(scorer) };

//     // 3) operators (simple defaults)
//     let mut selection = TournamentSelection::new(2, &mut rng);
//     let mut crossover  = DNFVecCrossover::new(&mut rng);
//     let mut mutation   = SimpleDNFVecMutation::new(/* go, gtex, conj_gen, rng */);
//     let elites         = ElitesByNumberSelector::new(1);

//     // 4) build GA (adjust to your actual constructor)
//     let mut ga = GeneticAlgorithm::new(
//         evaluator, Box::new(selection), Box::new(crossover),
//         Box::new(mutation), Box::new(elites),
//         /* phenotype: TermId, population_size, generations, mutation_rate, */ &mut rng,
//         /* formula_generator: */ Box::new(dnf_gen)
//     );

//     ga.initialize_population(4).unwrap();
//     let best = ga.fit();
//     assert!(best.get_score().is_finite());
// }
/* */
