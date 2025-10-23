use go2hpo_genetic_algorithm::annotations::GeneSetAnnotations;
use ontolius::ontology::OntologyTerms;
use ontolius::TermId;
use rand::{rng, Rng};
use rand::rngs::ThreadRng;

use rstest::rstest;

use go2hpo_genetic_algorithm::genetic_algorithm::{ConjunctionMutation, Mutation, Selection, SimpleDNFBitmaskMutation, SimpleDNFVecMutation};
use go2hpo_genetic_algorithm::genetic_algorithm::TournamentSelection;
use go2hpo_genetic_algorithm::genetic_algorithm::Crossover;
use go2hpo_genetic_algorithm::genetic_algorithm::{DNFVecCrossover, ConjunctionCrossover};
use go2hpo_genetic_algorithm::genetic_algorithm::ElitesSelector;
use go2hpo_genetic_algorithm::genetic_algorithm::ElitesByNumberSelector;

use go2hpo_genetic_algorithm::Solution;

use go2hpo_genetic_algorithm::logical_formula::{
    Conjunction, DNFBitmask, DNFVec, DgeState, Formula, RandomConjunctionGenerator, TermObservation, TissueExpression, DNF
};

use rand::{rngs::SmallRng, SeedableRng};
use gtex_analyzer::expression_analysis::GtexSummary;
use ontolius::ontology::csr::MinimalCsrOntology;
use go2hpo_genetic_algorithm::genetic_algorithm::GeneticAlgorithm;

use go2hpo_genetic_algorithm::utils::fixtures::gene_set_annotations::{
    go_sample,
    gtex_summary_sample,
    gene_set_annotations,
};


use go2hpo_genetic_algorithm::logical_formula::FormulaGenerator;
use go2hpo_genetic_algorithm::logical_formula::{GenePickerConjunctionGenerator, RandomDNFVecGenerator};



// ==========================
// 1) SELECTION
// ==========================
#[rstest]
fn selection_tournament_basic() {
    let mut rng: ThreadRng = rng();

    // dummy population of Solutions (just give them scores)
    let pop: Vec<Solution<usize>> = (0..5)
        .map(|i| Solution::new(i, i as f64)) // increasing score
        .collect();

    let mut sel = TournamentSelection::new(2, &mut rng);
    let picked = sel.select(&pop);
    // Must be one of the population members:
    assert!(pop.iter().any(|s| s.get_formula() == picked.get_formula()));
}

// ==========================
// 2) CROSSOVER (Conjunction)
// ==========================
#[rstest]
fn crossover_conjunction_balanced() {
    let mut rng: ThreadRng = rng();

    // Two tiny parents
    let p1 = Conjunction {
        term_observations: vec![
            TermObservation::new("GO:0051146".parse().unwrap(), true),
            TermObservation::new("GO:0052693".parse().unwrap(), false),
        ],
        tissue_expressions: vec![
            TissueExpression::new("Heart".to_string(), DgeState::Up),
        ],
    };

    let p2 = Conjunction {
        term_observations: vec![
            TermObservation::new("GO:0005634".parse().unwrap(), true),
        ],
        tissue_expressions: vec![
            TissueExpression::new("Liver".to_string(), DgeState::Down),
            TissueExpression::new("Lungs".to_string(), DgeState::Up),
        ],
    };

    // Balanced (parent1_fraction = 0.5 by default)
    let mut cx = ConjunctionCrossover::new(&mut rng);
    let child = cx.crossover(&p1, &p2);

    println!("{}", child);

    // Child should contain a subset mix of terms from both parents
    assert!(!child.term_observations.is_empty() || !child.tissue_expressions.is_empty());
}

// ==========================
// 3) CROSSOVER (DNFVec)
// ==========================
#[rstest]
fn crossover_dnfvec_balanced() {
    let mut rng: ThreadRng = rng();

    // Build two small DNFs (DNFVec holds concrete Conjunctions)
    let d1 = DNFVec::from_conjunctions(vec![
        Conjunction {
            term_observations: vec![TermObservation::new("GO:0005634".parse().unwrap(), true)],
            tissue_expressions: vec![],
        },
        Conjunction {
            term_observations: vec![TermObservation::new("GO:0051146".parse().unwrap(), false)],
            tissue_expressions: vec![],
        },
    ]);

    let d2 = DNFVec::from_conjunctions(vec![
        Conjunction {
            term_observations: vec![TermObservation::new("GO:0052693".parse().unwrap(), true), TermObservation::new("GO:0012345".parse().unwrap(), false)],
            tissue_expressions: vec![],
        },
        Conjunction {
            term_observations: vec![],
            tissue_expressions: vec![TissueExpression::new("Heart".to_string(), DgeState::Up)],
        },
        Conjunction {
            term_observations: vec![TermObservation::new("GO:0043215".parse().unwrap(), false)],
            tissue_expressions: vec![TissueExpression::new("Liver".to_string(), DgeState::Down)],
        }
    ]);

    let mut cx = DNFVecCrossover::new(&mut rng);
    let child = cx.crossover(&d1, &d2);
    println!("{}", child);

    // child should be a valid DNF with >= 1 conjunction from parents
    assert!(child.len() >= 1);
    // and its active conjunctions should be well-formed
    assert!(!child.get_active_conjunctions().is_empty());
}

// ==========================
// 4) SURVIVAL (elites only)
// ==========================
#[rstest]
fn survival_elites_by_number() {
    let mut next: Vec<Solution<usize>> = Vec::new();
    // Previous population with known scores
    let prev: Vec<Solution<usize>> = vec![
        Solution::new(10, 0.1),
        Solution::new(11, 0.9), // best
        Solution::new(12, 0.7),
        Solution::new(13, 0.3),
    ];

    let keep = 2;
    let surv = ElitesByNumberSelector::new(keep);
    let filled = surv.pass_elites(&mut next, &prev, false);

    assert_eq!(filled, keep);
    assert_eq!(next.len(), keep);
    // Ensure the best two by score are kept (0.9, 0.7)
    assert!(next.iter().any(|s| *s.get_formula() == 11));
    assert!(next.iter().any(|s| *s.get_formula() == 12));
}


#[rstest]
fn test_gene_picker_conjunction_generator(
    gene_set_annotations: GeneSetAnnotations,
) {
    
    let mut rng = SmallRng::seed_from_u64(42);
    let mut gen = GenePickerConjunctionGenerator::new(
            &mut rng,
            1.0,
            1.0,
            &gene_set_annotations,
            None,  // no HPO filter
            Some(3),
            Some(4)
        );

    // Generate a conjunction
    let conj: Conjunction = gen.generate();
    println!("Generated Conjunction: {}", conj);

    // Assertions: It should have at least one term or tissue expression
    // assert!(
    //     !conj.term_observations.is_empty() || !conj.tissue_expressions.is_empty(),
    //     "Conjunction should not be empty"
    // );

    // Ensure terms and tissues come from the gene set
    let annotations_map = gene_set_annotations.get_gene_annotations_map();
    let mut all_terms: Vec<_> = annotations_map.values().flat_map(|ga| ga.get_term_annotations()).collect();
    let mut all_tissues: Vec<_> = annotations_map.values().flat_map(|ga| ga.get_tissue_expressions()).collect();

    for obs in &conj.term_observations {
        assert!(all_terms.contains(&&obs.term_id), "Unexpected term {:?}", obs.term_id);
    }
    for tissue in &conj.tissue_expressions {
        assert!(all_tissues.contains(&tissue), "Unexpected tissue {:?}", tissue);
    }
}

#[rstest]
fn test_gene_picker_conjunction_generator_specific_hpo(
    gene_set_annotations: GeneSetAnnotations,
) {
    
    let mut rng = SmallRng::seed_from_u64(42);
    let target_hpo: TermId = "HP:0001250".parse().unwrap(); // e.g. Seizures
    let mut gen = GenePickerConjunctionGenerator::new(
        &mut rng,
        1.0,
        1.0,
        &gene_set_annotations,
        Some(target_hpo),  // restrict to this phenotype
        Some(2),
        Some(2)
    );

    // Generate a conjunction
    let conj: Conjunction = gen.generate();
    println!("Generated Conjunction: {}", conj);

    // Ensure terms and tissues come from the gene set
    let annotations_map = gene_set_annotations.get_gene_annotations_map();
    let mut all_terms: Vec<_> = annotations_map.values().flat_map(|ga| ga.get_term_annotations()).collect();
    let mut all_tissues: Vec<_> = annotations_map.values().flat_map(|ga| ga.get_tissue_expressions()).collect();

    for obs in &conj.term_observations {
        assert!(all_terms.contains(&&obs.term_id), "Unexpected term {:?}", obs.term_id);
    }
    for tissue in &conj.tissue_expressions {
        assert!(all_tissues.contains(&tissue), "Unexpected tissue {:?}", tissue);
    }
}




#[rstest]
fn test_random_dnfvec_generator(
    gene_set_annotations: GeneSetAnnotations,
) {
    let mut rng = SmallRng::seed_from_u64(99);
    let mut rng2 = SmallRng::seed_from_u64(99);

    // Create a simple GenePicker generator to feed into DNF generator
    let target_hpo: TermId = "HP:0001250".parse().unwrap(); // e.g. Seizures
    let mut conj_gen = GenePickerConjunctionGenerator::new(
        &mut rng,
        1.0,
        1.0,
        &gene_set_annotations,
        Some(target_hpo),  // restrict to this phenotype
        Some(2),
        Some(2)
    );

    // DNF generator with 3 conjunctions
    let mut dnf_gen = RandomDNFVecGenerator::new(&mut conj_gen, 3, rng2  );

    let dnf: DNFVec = dnf_gen.generate();
    println!("Generated DNFVec: {}", dnf);

    // Assertions
    let conjunctions = dnf.get_active_conjunctions();
    assert_eq!(conjunctions.len(), 3, "Expected exactly 3 conjunctions");

    // Ensure each conjunction has valid terms/tissues
    for conj in conjunctions {
        assert!(
            !conj.term_observations.is_empty() || !conj.tissue_expressions.is_empty(),
            "Each conjunction should not be empty"
        );
    }
}




#[rstest]
fn test_conjunction_mutation_sequence(
    go_sample: MinimalCsrOntology,
    gtex_summary_sample: std::io::Result<GtexSummary>,
) {
    let gtex = gtex_summary_sample.expect("Fixture GTEx summary must be ok");
    let mut rng = SmallRng::seed_from_u64(32);

    // Start with one GO term and one tissue expression
    let mut conj = Conjunction {
        term_observations: vec![TermObservation::new("GO:0051146".parse().unwrap(), false)],
        tissue_expressions: vec![TissueExpression::new("Liver".to_string(), DgeState::Up)],
    };

    let mut mutation = ConjunctionMutation::new(&go_sample, &gtex, 8, &mut rng);

    // 1) toggle_term_status
    println!("BEFORE TOGGLE TERM STATUS: {}", conj);
    let before_status = conj.term_observations[0].is_excluded;
    mutation.toggle_term_status(&mut conj);
    let after_status = conj.term_observations[0].is_excluded;
    println!("AFTER TOGGLE TERM STATUS: {}", conj);
    assert_ne!(before_status, after_status);

    // 2) add_random_term
    let before_len_terms = conj.term_observations.len();
    mutation.add_random_term(&mut conj);
    let after_len_terms = conj.term_observations.len();
    println!("AFTER ADD RANDOM TERM: {}", conj);
    assert_eq!(after_len_terms, before_len_terms + 1);

    // 3) delete_random_term
    let before_len_terms = conj.term_observations.len();
    mutation.delete_random_term(&mut conj);
    let after_len_terms = conj.term_observations.len();
    println!("AFTER DELETE RANDOM TERM: {}", conj);
    assert_eq!(after_len_terms, before_len_terms.saturating_sub(1));

    // 4) toggle_tissue_expression_state
    let before_state = conj.tissue_expressions[0].state.clone();
    mutation.toggle_tissue_expression_state(&mut conj);
    let after_state = conj.tissue_expressions[0].state.clone();
    println!("AFTER TOGGLE TISSUE EXPRESSION STATE: {}", conj);
    assert_ne!(before_state, after_state);

    // 5) add_tissue_expression_term
    let before_len_tissues = conj.tissue_expressions.len();
    mutation.add_tissue_expression_term(&mut conj);
    let after_len_tissues = conj.tissue_expressions.len();
    println!("AFTER ADD TISSUE EXPRESSION TERM: {}", conj);
    assert_eq!(after_len_tissues, before_len_tissues + 1);

    // 6) delete_tissue_expression_term
    let before_len_tissues = conj.tissue_expressions.len();
    mutation.delete_tissue_expression_term(&mut conj);
    let after_len_tissues = conj.tissue_expressions.len();
    println!("AFTER DELETE TISSUE EXPRESSION TERM: {}", conj);
    assert_eq!(after_len_tissues, before_len_tissues.saturating_sub(1));
}


#[rstest]
fn test_conjunction_mutation_parent_child(
    go_sample: MinimalCsrOntology,
    gtex_summary_sample: std::io::Result<GtexSummary>,
) {
    let gtex = gtex_summary_sample.expect("Fixture GTEx summary must be ok");
    let mut rng = SmallRng::seed_from_u64(42);

    // Take a few GO and tissue terms
    let go_terms: Vec<_> = go_sample.iter_term_ids().take(5).cloned().collect();
    let tissue_terms: Vec<String> = gtex.metadata.get_tissue_names().into_iter().cloned().take(5).collect();

    // Create a conjunction with one GO term
    let mut conj = Conjunction {
        term_observations: vec![
            TermObservation::new(go_terms[0].clone(), false),
        ],
        tissue_expressions: vec![],
    };

    let mut conj2 = Conjunction {
        term_observations: vec![
            TermObservation::new(go_terms[2].clone(), false),
        ],
        tissue_expressions: vec![],
    };

    let mut mutation = ConjunctionMutation::new(&go_sample, &gtex, 8, &mut rng);

    // CONJUNCTION 1
    println!("INITIAL CONJUNCTION: {}", conj);
    // Test mutate_with_parent_term 
    mutation.mutate_with_parent_term(&mut conj);
    println!("AFTER PARENT MUTATION: {}", conj);
    assert!(conj.term_observations.len() >= 1);

    // Test mutate_with_child_term
    mutation.mutate_with_child_term(&mut conj);
    println!("AFTER CHILD MUTATION: {}", conj);
    assert!(conj.term_observations.len() >= 1);

    for obs in &conj.term_observations {
        let exists = go_sample.iter_term_ids().any(|tid| tid == &obs.term_id);
        assert!(exists, "Term {:?} is not in ontology", obs.term_id);
    }


    // CONJUNCTION 2
    println!("\nINITIAL CONJUNCTION: {}", conj2);
    // Test mutate_with_parent_term 
    mutation.mutate_with_parent_term(&mut conj2);
    println!("AFTER PARENT MUTATION: {}", conj2);
    assert!(conj2.term_observations.len() >= 1);

    // Test mutate_with_child_term
    mutation.mutate_with_child_term(&mut conj2);
    println!("AFTER CHILD MUTATION: {}", conj2);
    assert!(conj2.term_observations.len() >= 1);

    for obs in &conj2.term_observations {
        let exists = go_sample.iter_term_ids().any(|tid| tid == &obs.term_id);
        assert!(exists, "Term {:?} is not in ontology", obs.term_id);
    }


}




#[rstest]
fn test_dnfbitmask_mutation_toggle(
    go_sample: MinimalCsrOntology,
    gtex_summary_sample: std::io::Result<GtexSummary>
) {
    let gtex = gtex_summary_sample.expect("Fixture GTEx summary must be ok");
    let mut rng = SmallRng::seed_from_u64(7);

    // build a trivial DNFBitmask with 2 conjunctions
    let conj1 = Conjunction { term_observations: vec![TermObservation::new("GO:0051146".parse().unwrap(), false)], tissue_expressions: vec![] };
    let conj2 = Conjunction { term_observations: vec![], tissue_expressions: vec![TissueExpression::new("Liver".to_string(), DgeState::Up)] };
    let conj3 = Conjunction { term_observations: vec![TermObservation::new("GO:0012345".parse().unwrap(), true)], tissue_expressions: vec![TissueExpression::new("Heart".to_string(), DgeState::Up)] };
    let conjunctions = [conj1, conj2, conj3];
    let mut dnf = DNFBitmask::new(&conjunctions);
    dnf.activate_conjunction(1).expect("It should activate the conjunction at position 1");
    println!("BEFORE BITMASK MUTATION: {}",dnf);
    let before = dnf.len();
    let mut mutator = SimpleDNFBitmaskMutation::new(&mut rng);
    mutator.mutate(&mut dnf);
    println!("AFTER BITMASK FIRST MUTATION: {}",dnf);
    let after = dnf.len();
    assert!(after == before + 1 || after == before.saturating_sub(1));

    mutator.mutate(&mut dnf);
    println!("AFTER BITMASK SECOND MUTATION {}",dnf);
    let after2 = dnf.len();
    assert!(after2 == after + 1 || after2 == after.saturating_sub(1));
}


#[rstest]
fn test_simple_dnfvec_mutation(
    go_sample: MinimalCsrOntology,
    gtex_summary_sample: std::io::Result<GtexSummary>,
) {
    let gtex = gtex_summary_sample.expect("Fixture GTEx summary must be ok");
    let mut rng = SmallRng::seed_from_u64(123);
    let mut rng2 = SmallRng::seed_from_u64(30); 

    // Build a random Conjunction generator using a small GO + GTEx subset
    let go_terms: Vec<_> = go_sample.iter_term_ids().take(5).cloned().collect();
    let tissue_terms: Vec<String> = gtex.metadata.get_tissue_names().into_iter().cloned().take(5).collect();

    let conj_gen = RandomConjunctionGenerator::new(2, &go_terms, 2, &tissue_terms, rng.clone());

    // Create the ConjunctionMutation and SimpleDNFVecMutation
    let conj_mut = ConjunctionMutation::new(&go_sample, &gtex, 8, &mut rng);
    let mut dnf_mut = SimpleDNFVecMutation::new(conj_mut, conj_gen, 4, &mut rng2);
    // Initial DNFVec with 2 conjunctions
    let mut dnf = DNFVec::from_conjunctions(vec![
        Conjunction {
            term_observations: vec![TermObservation::new("GO:0051146".parse().unwrap(), false)],
            tissue_expressions: vec![],
        },
         Conjunction { 
            term_observations: vec![TermObservation::new("GO:0005634".parse().unwrap(), true)], 
            tissue_expressions: vec![TissueExpression::new("Heart".to_string(), DgeState::Up)] 
        },
    ]);
    println!("BEFORE VECTOR MUTATION: {}",dnf);
    let before_len = dnf.len();
    dnf_mut.mutate(&mut dnf);
    let after_len = dnf.len();
    println!("AFTER FIRST VECTOR MUTATION: {}",dnf);
    dnf_mut.mutate(&mut dnf);
    let after_len2 = dnf.len();
    println!("AFTER SECOND VECTOR MUTATION: {}",dnf);
    dnf_mut.mutate(&mut dnf);
    let after_len3 = dnf.len();
    println!("AFTER THIRD VECTOR MUTATION: {}",dnf);

}

use go2hpo_genetic_algorithm::logical_formula::NaiveSatisfactionChecker;
use go2hpo_genetic_algorithm::genetic_algorithm::{ConjunctionScorer, DNFScorer, FormulaEvaluator, ScoreMetric};

#[rstest]
fn test_formula_evaluator_with_conjunctionscorer(
    go_sample: MinimalCsrOntology,
    gene_set_annotations: GeneSetAnnotations,
) {


    let checker = NaiveSatisfactionChecker::new(&go_sample, &gene_set_annotations);
    let scorer = ConjunctionScorer::new(checker, ScoreMetric::Accuracy);
    let evaluator = FormulaEvaluator::new(Box::new(scorer));

    let t1: TermId = "GO:0051146".parse().unwrap();
    let conj = Conjunction {
        term_observations: vec![TermObservation::new(t1.clone(), false)],
        tissue_expressions: vec![TissueExpression::new("Liver".to_string(), DgeState::Up)],
    };

    println!("{}", conj);
    let phenotype: TermId = "HPO:0000001".parse().unwrap();

    let solution = evaluator.evaluate(&conj, &phenotype);
    println!("Conjunction Solution: {}", solution);

    assert!(solution.get_score() >= 0.0 && solution.get_score() <= 1.0);
}

#[rstest]
fn test_formula_evaluator_with_dnfscorer(
    go_sample: MinimalCsrOntology,
    gene_set_annotations: GeneSetAnnotations,
) {


    let checker = NaiveSatisfactionChecker::new(&go_sample, &gene_set_annotations);
    let conjunction_scorer = ConjunctionScorer::new(checker, ScoreMetric::Accuracy);
    let scorer = DNFScorer::new(conjunction_scorer, 0.3);
    let evaluator = FormulaEvaluator::new(Box::new(scorer));;

    let t1: TermId = "GO:0051146".parse().unwrap();
    let conj = Conjunction {
        term_observations: vec![TermObservation::new(t1.clone(), false)],
        tissue_expressions: vec![TissueExpression::new("Brain".to_string(), DgeState::Down)],
    };

    let dnf = DNFVec::from_conjunctions(vec![conj]);
    println!("{}", dnf);
    let phenotype: TermId = "HPO:0000001".parse().unwrap();

    let solution = evaluator.evaluate(&dnf, &phenotype);
    println!("DNF Solution: {}", solution);

    assert!(solution.get_score() >= 0.0 && solution.get_score() <= 1.0);
}



#[rstest]
fn test_genetic_algorithm_sanity(
    go_sample: MinimalCsrOntology,
    gene_set_annotations: GeneSetAnnotations,
    gtex_summary_sample: std::io::Result<GtexSummary>,
) {
    let gtex = gtex_summary_sample.expect("Fixture GTEx summary must be ok");

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
    let checker = NaiveSatisfactionChecker::new(&go_sample, &gene_set_annotations);
    let conj_scorer = ConjunctionScorer::new(checker, ScoreMetric::Accuracy);
    let scorer = DNFScorer::new(conj_scorer, 0.0);
    let evaluator = FormulaEvaluator::new(Box::new(scorer));

    // Operators (keep simple)
    let go_terms: Vec<_> = go_sample.iter_term_ids().take(5).cloned().collect();
    let tissue_terms: Vec<String> = gtex.metadata.get_tissue_names().into_iter().cloned().take(5).collect();
    let selection = Box::new(TournamentSelection::new(2, &mut rng_selection));
    let crossover = Box::new(DNFVecCrossover::new(&mut rng_crossover));
    let mutation = Box::new(SimpleDNFVecMutation::new(
        ConjunctionMutation::new(&go_sample, &gtex, 8, &mut rng_conj_mutation),
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

    // Build GA with just 1 generation
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
        "HP:0000001".parse().unwrap(),  // dummy phenotype
    );

    // Run one generation
    let best = ga.fit();
    println!("Sanity GA best solution: {}", best);

    // Assertions: population size stays the same, score is finite
    assert_eq!(ga.get_population().len(), pop_size);
    assert!(best.get_score().is_finite());
}





// #[rstest]
// fn test_genetic_algorithm_end_to_end(
//     go_sample: MinimalCsrOntology,
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

//     let terms_mutation_pool: Vec<TermId> = go_sample.iter_term_ids().take(10).cloned().collect();
//     let tissues_mutation_pool: Vec<String> = gtex.metadata.get_tissue_names().into_iter().take(10).cloned().collect();

//     let mut dnf_gen = RandomDNFVecGenerator::new(&mut conj_gen, 3, rng.clone());

//     // 2) Build evaluator (ConjunctionScorer -> DNFScorer -> FormulaEvaluator)
//     let checker = NaiveSatisfactionChecker::new(go_sample, &gene_set_annotations);
//     let conj_scorer = ConjunctionScorer::new(checker, ScoreMetric::Accuracy);
//     let scorer = DNFScorer::new(conj_scorer);
//     let evaluator = FormulaEvaluator::new(Box::new(scorer));

//     // 3) Operators
//     let selection = Box::new(TournamentSelection::new(2, &mut rng2));
//     let crossover = Box::new(DNFVecCrossover::new(&mut rng3));
//     let mutation = Box::new(SimpleDNFVecMutation::new(
//         ConjunctionMutation::new(&go_sample, &gtex, 8, &mut rng),
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
