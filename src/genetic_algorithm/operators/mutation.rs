use crate::genetic_algorithm::Solution;
use rand::prelude::*;
use crate::logical_formula::{ConjunctionGenerator, DgeState, Formula, TissueExpression};
use gtex_analyzer::expression_analysis::GtexSummary;

use ontolius::term::simple::SimpleMinimalTerm;
use ontolius::{
    ontology::{HierarchyWalks, OntologyTerms},
    TermId,
};


use crate::logical_formula::DNFVec;
use crate::{
    logical_formula::{Conjunction, TermObservation},
    logical_formula::{DNFBitmask, DNF},
};


// FORMULA MUTATION OPERATOR
pub trait Mutation<T> {
    fn mutate(&mut self, formula: &mut T);
}

pub struct ConjunctionMutation<'a, O, R: Rng> {
    go: &'a O,
    gtex: &'a GtexSummary,
    max_n_terms: usize,
    rng: &'a mut R,
}

impl<O, R> Mutation<Conjunction> for ConjunctionMutation<'_, O, R>
where
    O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
    R: Rng,
{
    fn mutate(&mut self, formula: &mut Conjunction) {
        let rnd_num = self.rng.random_range(0..=7);
        match rnd_num {
            0 => self.mutate_with_parent_term(formula),
            1 => self.mutate_with_child_term(formula),
            2 => self.delete_random_term(formula),
            3 => self.add_random_term(formula),
            4 => self.toggle_term_status(formula),
            5 => self.delete_tissue_expression_term(formula),
            6 => self.add_tissue_expression_term(formula),
            7 => self.toggle_tissue_expression_state(formula),
            _ => panic!("A random number outside of the range has been generated. No associated mutation"),
        }
    }
}

impl<'a, O, R> ConjunctionMutation<'a, O, R>
where
    O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
    R: Rng,
{
    /// Exchange a GO term with a parent
    pub fn mutate_with_parent_term(&mut self, formula: &mut Conjunction) {
        // Get a term from a random index. Maybe a more sophisticated way will be used in the future
        if formula.term_observations.is_empty() {
            self.add_random_term(formula);
            return; 
        }
        let rnd_index = self.rng.random_range(0..formula.term_observations.len());
        let mut term_ob = formula
            .term_observations
            .get_mut(rnd_index)
            .expect("It should return a term");

        // Get the parents' term IDs, I collect them to know the length without consuming it
        let parents_ids: Vec<&TermId> = self.go.iter_parent_ids(&term_ob.term_id).collect(); //or go.iter_ancestor_ids(term)

        if parents_ids.len() > 0{
            // Select one of them randomly and subsitute it with the current termId in the formula
            let rnd_index = self.rng.random_range(0..parents_ids.len());
            term_ob.term_id = parents_ids.get(rnd_index).copied().unwrap().clone();
        }
    }

    ///  Exchange a GO term with one of its children
    pub fn mutate_with_child_term(&mut self, formula: &mut Conjunction) {
        // Get a term from a random index. Maybe a more sophisticated way will be used in the future
        if formula.term_observations.is_empty() {
            return;
        }
        let rnd_index = self.rng.random_range(0..formula.term_observations.len());
        let mut term_ob = formula
            .term_observations
            .get_mut(rnd_index)
            .expect("It should return a term");

        // Get the childrens' term IDs, I collect them to know the length without consuming it
        let child_ids: Vec<&TermId> = self.go.iter_child_ids(&term_ob.term_id).collect(); //or go.iter_descendant_ids(term)

        if child_ids.len() > 0{   
            // Select one of them randomly and subsitute it with the current termId in the formula
            let rnd_index = self.rng.random_range(0..child_ids.len());
            term_ob.term_id = child_ids.get(rnd_index).copied().unwrap().clone();
        }
    }

    /// Delete a GO term
    pub fn delete_random_term(&mut self, formula: &mut Conjunction) {
        if formula.term_observations.is_empty() {
            return; // or handle with Result
        }
        let rnd_index = self.rng.random_range(0..formula.term_observations.len());
        formula.term_observations.remove(rnd_index);
    }

    /// Add a random GO term
    pub fn add_random_term(&mut self, formula: &mut Conjunction) {
        if formula.len() >= self.max_n_terms{
            let rnd_indx = self.rng.random_range(0..4);
            match rnd_indx {
                0 => self.mutate_with_parent_term(formula),
                1 => self.mutate_with_child_term(formula),
                2 => self.delete_random_term(formula),
                3 => self.toggle_term_status(formula),
                _ => panic!("A random number outside of the range has been generated in add_random_term"),
            };
            return;
        }
        let rnd_index = self.rng.random_range(0..self.go.len());
        if let Some(new_term) = self.go.iter_term_ids().nth(rnd_index) {
            let term_obs = TermObservation::new(new_term.clone(), self.rng.random_bool(0.5));
            formula.term_observations.push(term_obs);
        }
    }

    /// Toggle the status of GO terms with is_excluded
    pub fn toggle_term_status(&mut self, formula: &mut Conjunction) {
        if formula.term_observations.is_empty() {
            return;
        }
        let rnd_index = self.rng.random_range(0..formula.term_observations.len());
        let mut term_ob = formula
            .term_observations
            .get_mut(rnd_index)
            .expect("It should return a term");
        term_ob.is_excluded = !term_ob.is_excluded;
    }

    /// Delete a gene expression term
    pub fn delete_tissue_expression_term(&mut self, formula: &mut Conjunction) {
        
        if formula.tissue_expressions.is_empty() {
            return; // or handle with Result
        }
        let rnd_index = self.rng.random_range(0..formula.tissue_expressions.len());
        formula.tissue_expressions.remove(rnd_index);
    }

    /// Add a gene expression term from a random tissue
    pub fn add_tissue_expression_term(&mut self, formula: &mut Conjunction) {
        if formula.len() >= self.max_n_terms{
            let rnd_indx = self.rng.random_range(0..2);
            match rnd_indx {
                0 => self.toggle_tissue_expression_state(formula),
                1 => self.delete_tissue_expression_term(formula),
                _ => panic!("A random number outside of the range has been generated in add_tissue_expression_term"),
            };
            return;
        }

        let tissues = self.gtex.metadata.get_tissue_names();
        
        let rnd_index = self.rng.random_range(0..tissues.len());
        
        if let Some(tissue_name) = tissues.get(rnd_index){
            let tissue_expr = TissueExpression::new(tissue_name.clone(), DgeState::get_random(self.rng));
            formula.tissue_expressions.push(tissue_expr)
        }
    }

    // toggle of a tissue expression state from lo to hi or vice versa
    pub fn toggle_tissue_expression_state(&mut self, formula: &mut Conjunction) {
        if formula.tissue_expressions.is_empty() {
            return; // or handle with Result
        }
        let rnd_index = self.rng.random_range(0..formula.tissue_expressions.len());
        let mut tissue_term = formula.tissue_expressions.get_mut(rnd_index).expect("It should return a TissueExpression");
        
        let dge_states = [DgeState::Down, DgeState::Normal, DgeState::Up];

        let new_state = match tissue_term.state {
            DgeState::Down =>  [DgeState::Normal, DgeState::Up].choose(&mut self.rng).expect("Should return one state"),
            DgeState::Normal =>  [DgeState::Down, DgeState::Up].choose(&mut self.rng).expect("Should return one state"),
            DgeState::Up =>  [DgeState::Normal, DgeState::Down].choose(&mut self.rng).expect("Should return one state"),
        };

        tissue_term.state =  new_state.clone();
    }

    
    pub fn new(go: &'a O, gtex: &'a GtexSummary, max_n_terms: usize, rng: &'a mut R,) -> Self{
        Self{go, gtex, max_n_terms, rng}
    }
    
}

pub struct SimpleDNFBitmaskMutation<'a, R: Rng>{
    rng: &'a mut R,
}

impl<'a, R: Rng> Mutation<DNFBitmask<'_>> for SimpleDNFBitmaskMutation<'a, R> {
    fn mutate(&mut self, formula: &mut DNFBitmask) {
        let n_conjunctions = formula.total_conjunctions_count();
        let rnd_index = self.rng.random_range(0..n_conjunctions);

        formula
            .toggle_conjunction(rnd_index)
            .expect("It should toggle the conjunction");
    }
}

impl<'a, R: Rng> SimpleDNFBitmaskMutation<'a, R> {
    pub fn new(rng: &'a mut R) -> Self{
        Self {rng}
    }
}


pub struct SimpleDNFVecMutation<'a, O, G, R> 
where 
O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
G: ConjunctionGenerator,
R: Rng,
{
    conjunction_mutation: ConjunctionMutation<'a, O, R>,
    conjunction_generator: G,
    max_n_conj: usize,
    rng: &'a mut R,
}

impl<'a, O, G, R> SimpleDNFVecMutation<'a, O, G, R> 
where 
O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
G: ConjunctionGenerator,
R: Rng,{
    pub fn new(
        conjunction_mutation: ConjunctionMutation<'a, O, R>,
        conjunction_generator: G,
        max_n_conj: usize,
        rng: &'a mut R,
    ) -> Self {
        Self {
            conjunction_mutation,
            conjunction_generator,
            max_n_conj,
            rng,
        }
    }
}


impl<'a, O, G, R> Mutation<DNFVec> for SimpleDNFVecMutation<'a, O, G, R>
where
    O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
    G: ConjunctionGenerator,
    R: Rng, {
    fn mutate(&mut self, formula: &mut DNFVec) {
        let rnd_num = self.rng.random_range(0..=2);
        match rnd_num {
            0 => self.mutate_conjunction(formula),
            1 => self.add_random_conjunction(formula),
            2 => self.remove_random_conjunction(formula),
            _ => panic!("A random number outside of the range has been generated. No associated mutation"),
        }
    }
}


impl<'a, O, G, R> SimpleDNFVecMutation<'a, O, G, R>
where
    O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
    G: ConjunctionGenerator,
    R: Rng,
{
    /// Choose one conjunction randomly and then modifies it with ConjunctionMutation
    pub fn mutate_conjunction(&mut self, formula: &mut DNFVec) {
        let mut conjunctions = formula.get_mut_active_conjunctions();
        let rnd_index = self.rng.random_range(0..conjunctions.len());
        let mut conjunction = conjunctions
            .get_mut(rnd_index)
            .expect("It should return a conjunction");
        self.conjunction_mutation.mutate(&mut conjunction);
    }

    /// creates a new conjunction and assigns it to the formula
    pub fn add_random_conjunction(&mut self, formula: &mut DNFVec) {
        // If the max number of conjunctions has been reached, instead of adding a conjunction, mutate 1
        if formula.len() >= self.max_n_conj{
            self.mutate_conjunction(formula);
            return;
        }
        let conjunction: Conjunction = self.conjunction_generator.generate();
        formula.activate_conjunction(conjunction).expect("The conjunction should be added without errors")
    }

    /// Removes a conjunction randomly 
    pub fn remove_random_conjunction(&mut self, formula: &mut DNFVec) {
        let conjunctions = formula.get_mut_active_conjunctions();

        if conjunctions.is_empty() {
            self.add_random_conjunction(formula);
            return;
        }

        let rnd_index = self.rng.random_range(0..conjunctions.len());
        conjunctions.remove(rnd_index);
    }

}

pub struct BiasedDNFMutation;

impl Mutation<DNFBitmask<'_>> for BiasedDNFMutation {
    fn mutate(&mut self, formula: &mut DNFBitmask) {
        todo!()
    }
}


#[cfg(test)]
mod tests {
    use lazy_static::lazy_static;
    use num::Saturating;

    use crate::logical_formula::{Formula, FormulaGenerator, RandomConjunctionGenerator, RandomDNFBistmaskGenerator};

    use super::*;

    use std::{fs::File, io::BufReader};

    use flate2::bufread::GzDecoder;
    use ontolius::{
        io::OntologyLoaderBuilder,
        ontology::{csr::{CsrOntology, MinimalCsrOntology}, HierarchyQueries},
    };

    // use ontolius::ontology::OntologyTerms;
    use gtex_analyzer::expression_analysis::GtexSummaryLoader;


    lazy_static! {
        // static ref GO: ontolius::ontology::csr::CsrOntology<u32, SimpleMinimalTerm> = get_go(); 
        static ref small_test_tissues: Vec<String> = vec!["Colon_Transverse_Muscularis".to_string(),
                                                         "Colon_Transverse_Mixed_Cell".to_string(),
                                                          "Colon_Transverse_Muscularis".to_string(),
                                                          "Testis".to_string(),
                                                          "Small_Intestine_Terminal_Ileum_Mixed_Cell".to_string()];
        static ref small_test_terms: Vec<TermId> = vec!["GO:0048856","GO:0110165","GO:0051146", "GO:0042692","GO:0052693", "GO:0030154","GO:0005634"]
            .into_iter()
            .map(|s| s.parse().unwrap())
            .collect();
        static ref GTEX: GtexSummary =  get_gtex_summary();
        static ref SMALL_TISS_EXPR_LIST: Vec<TissueExpression> = get_tissue_express_vec(30, 10);
        static ref SMALL_TERM_OBS_LIST: Vec<TermObservation> = get_term_obs_vec(30, 10, small_test_terms.to_vec());
        static ref SEED_LIST: Vec<u64> = vec![9, 12, 15, 18, 30];
        static ref SMALL_CONJUNCTION_LIST: Vec<Conjunction> = get_random_conjunction_list(10, 2, &small_test_terms.to_vec(), 2, &small_test_tissues.to_vec(), 64);
    }

    fn get_go() -> MinimalCsrOntology{
        let go_path = "data/go/go.toy.json.gz";
        let reader = GzDecoder::new(BufReader::new(
            File::open(go_path).expect("The file should be in the repo"),
        ));

        let parser = OntologyLoaderBuilder::new().obographs_parser().build();
        let go: MinimalCsrOntology = parser
            .load_from_read(reader)
            .expect("The ontology file should be OK");
        go
    }

    fn get_gtex_summary() -> GtexSummary{
        let file_path: &str = "data/gtex/GTEx_RNASeq_gene_median_tpm_HEAD.gct";

        let file = File::open(file_path).expect("The gtex file should open");
        let reader = BufReader::new(file);

        let summary_loader = GtexSummaryLoader::new(Some(10), None);
        let summary = summary_loader.load_summary(reader);
        summary.unwrap()
    }

    fn get_tissue_express_vec(seed: u64, size: usize) -> Vec<TissueExpression>{
        let mut uber_rng = SmallRng::seed_from_u64(seed);
        let tiss_exprs : Vec<TissueExpression> = (0..size)
                        .map(|i| (i,DgeState::get_random(&mut uber_rng)))
                        .map(|(i, state)| TissueExpression::new(format!("tissue_{}", i), state))
                        .collect();
        tiss_exprs
    }

    fn get_term_obs_vec(seed: u64, size: usize, term_ids: Vec<TermId>) -> Vec<TermObservation>{
        let mut uber_rng = SmallRng::seed_from_u64(seed);
        let mut term_chooser = SmallRng::seed_from_u64(seed);
        let term_obs : Vec<TermObservation> = (0..size)
                        .map(|i| (i, uber_rng.random_bool(0.5)))
                        .map(|(i, state)| TermObservation::new(term_ids.choose(&mut term_chooser).unwrap().clone(), state))
                        .collect();
        term_obs
    }

    fn get_random_conjunction_list<'a>(list_size: usize,n_go_terms: usize, go_terms: &'a Vec<TermId>, n_tissue_terms: usize, tissue_terms: &'a Vec<String>, seed: u64) -> Vec<Conjunction>{
        let mut uber_rng = SmallRng::seed_from_u64(seed);
        let mut random_conjunction_generator = RandomConjunctionGenerator::new(n_go_terms, go_terms, n_tissue_terms, tissue_terms, uber_rng);
        let small_conjunction_list: Vec<Conjunction> = (0..list_size).map(|_| random_conjunction_generator.generate()).collect(); 
        small_conjunction_list
    }

    #[test]
    fn test_toggle_tissue_expression_state(){
        let go = get_go(); 
        for seed in SEED_LIST.to_vec(){
            let mut formula = Conjunction::new();
            formula.tissue_expressions = SMALL_TISS_EXPR_LIST.clone();
            let mut rng = SmallRng::seed_from_u64(seed);
            let mut conjunction_mutation = ConjunctionMutation::new(&go, &GTEX, 8, &mut rng);
    
            let mut rng_twin = SmallRng::seed_from_u64(seed);
            let rnd_index = rng_twin.random_range(0..formula.tissue_expressions.len());
            
            conjunction_mutation.toggle_tissue_expression_state(&mut formula);
            
            for i in (0..SMALL_TISS_EXPR_LIST.len()){
                if i == rnd_index{
                    assert_ne!(SMALL_TISS_EXPR_LIST.get(rnd_index).expect("It should be Some").state,formula.tissue_expressions.get(rnd_index).expect("it should be Some").state)
                }else{
                    assert_eq!(SMALL_TISS_EXPR_LIST.get(i).expect("It should be Some").state,formula.tissue_expressions.get(i).expect("it should be Some").state)
                }
            }
        }
        
    }

    #[test]
    fn test_add_tissue_expression_term(){
        let go = get_go(); 
        for seed in SEED_LIST.to_vec(){
            let mut formula = Conjunction::new();
            formula.tissue_expressions = SMALL_TISS_EXPR_LIST.clone();
            let mut rng = SmallRng::seed_from_u64(seed);
            let mut conjunction_mutation = ConjunctionMutation::new(&go, &GTEX, 6, &mut rng);
            
            conjunction_mutation.add_tissue_expression_term(&mut formula);
            
            assert_eq!(SMALL_TISS_EXPR_LIST.len() + 1, formula.tissue_expressions.len())
        }
    }

    #[test]
    fn test_delete_tissue_expression_term(){
        let go = get_go(); 
        for seed in SEED_LIST.to_vec(){
            let mut formula = Conjunction::new();
            formula.tissue_expressions = SMALL_TISS_EXPR_LIST.clone();
            let mut rng = SmallRng::seed_from_u64(seed);
            let mut conjunction_mutation = ConjunctionMutation::new(&go, &GTEX, 10, &mut rng);
            
            conjunction_mutation.delete_tissue_expression_term(&mut formula);
            
            assert_eq!(SMALL_TISS_EXPR_LIST.len().saturating_sub(1), formula.tissue_expressions.len())
        }
    }

    #[test]
    fn test_toggle_term_status(){
        let go: MinimalCsrOntology = get_go(); 
        for seed in SEED_LIST.to_vec(){
            let mut formula = Conjunction::new();
            formula.term_observations = SMALL_TERM_OBS_LIST.clone();
            let mut rng = SmallRng::seed_from_u64(seed);
            let mut conjunction_mutation = ConjunctionMutation::new(&go, &GTEX, 12, &mut rng);
    
            let mut rng_twin = SmallRng::seed_from_u64(seed);
            let rnd_index = rng_twin.random_range(0..formula.term_observations.len());
            
            conjunction_mutation.toggle_term_status(&mut formula);
            
            for i in (0..SMALL_TERM_OBS_LIST.len()){
                if i == rnd_index{
                    assert_ne!(SMALL_TERM_OBS_LIST.get(rnd_index).expect("It should be Some").is_excluded, formula.term_observations.get(rnd_index).expect("it should be Some").is_excluded)
                }else{
                    assert_eq!(SMALL_TERM_OBS_LIST.get(i).expect("It should be Some").is_excluded, formula.term_observations.get(i).expect("it should be Some").is_excluded)
                }
            }
        }
    }

    #[test]
    fn test_add_random_term(){
        let go = get_go(); 
        for seed in SEED_LIST.to_vec(){
            let mut formula = Conjunction::new();
            formula.term_observations = SMALL_TERM_OBS_LIST.clone();
            let mut rng = SmallRng::seed_from_u64(seed);
            let mut conjunction_mutation = ConjunctionMutation::new(&go, &GTEX, 8, &mut rng);
                
            conjunction_mutation.add_random_term(&mut formula);
            
            assert_eq!(SMALL_TERM_OBS_LIST.len() + 1, formula.term_observations.len())
        }
    }

    #[test]
    fn test_delete_random_term(){
        let go = get_go(); 
        for seed in SEED_LIST.to_vec(){
            let mut formula = Conjunction::new();
            formula.term_observations = SMALL_TERM_OBS_LIST.clone();
            let mut rng = SmallRng::seed_from_u64(seed);
            let mut conjunction_mutation = ConjunctionMutation::new(&go, &GTEX, 10, &mut rng);
            
            conjunction_mutation.delete_random_term(&mut formula);
            
            assert_eq!(SMALL_TERM_OBS_LIST.len().saturating_sub(1), formula.term_observations.len())
        }
    }

    #[test]
    fn test_mutate_with_child_term(){
        let go: MinimalCsrOntology = get_go(); 
        for seed in SEED_LIST.to_vec(){
            let mut formula = Conjunction::new();
            formula.term_observations = SMALL_TERM_OBS_LIST.clone();
            let mut rng = SmallRng::seed_from_u64(seed);
            let mut conjunction_mutation = ConjunctionMutation::new(&go, &GTEX, 10, &mut rng);
    
            let mut rng_twin = SmallRng::seed_from_u64(seed);
            let rnd_index = rng_twin.random_range(0..formula.term_observations.len());
            
            conjunction_mutation.mutate_with_child_term(&mut formula);
            
            for i in (0..SMALL_TERM_OBS_LIST.len()){
                if i == rnd_index{
                    let parent = SMALL_TERM_OBS_LIST.get(rnd_index).expect("It should be Some").term_id.clone();
                    let children_ids: Vec<&TermId> = go.iter_child_ids(&parent).collect();
                    if children_ids.len() > 0{
                        let child = formula.term_observations.get(rnd_index).expect("It should be Some").term_id.clone();
                        assert!(go.is_child_of(&child, &parent));
                    }
                }else{
                    assert_eq!(SMALL_TERM_OBS_LIST.get(i).expect("It should be Some"), formula.term_observations.get(i).expect("it should be Some"))
                }
            }
        }
    }

    #[test]
    fn test_mutate_with_parent_term(){
        let go: MinimalCsrOntology = get_go(); 
        for seed in SEED_LIST.to_vec(){
            let mut formula = Conjunction::new();
            formula.term_observations = SMALL_TERM_OBS_LIST.clone();
            let mut rng = SmallRng::seed_from_u64(seed);
            let mut conjunction_mutation = ConjunctionMutation::new(&go, &GTEX, 10, &mut rng);
    
            let mut rng_twin = SmallRng::seed_from_u64(seed);
            let rnd_index = rng_twin.random_range(0..formula.term_observations.len());
            
            conjunction_mutation.mutate_with_parent_term(&mut formula);
            
            for i in (0..SMALL_TERM_OBS_LIST.len()){
                if i == rnd_index{
                    let child = SMALL_TERM_OBS_LIST.get(rnd_index).expect("It should be Some").term_id.clone();
                    let parent_ids: Vec<&TermId> = go.iter_parent_ids(&child).collect();
                    if parent_ids.len() > 0{
                        let parent = formula.term_observations.get(rnd_index).expect("It should be Some").term_id.clone();
                        assert!(go.is_parent_of(&parent, &child));
                    }
                }else{
                    assert_eq!(SMALL_TERM_OBS_LIST.get(i).expect("It should be Some"), formula.term_observations.get(i).expect("it should be Some"))
                }
            }
        }
    }

    #[test]
    fn test_mutate_dnfbitmask(){
        let go: MinimalCsrOntology = get_go();
        // let conjunction_generator = 
        let conjunctions = SMALL_CONJUNCTION_LIST.to_vec();
        for seed in SEED_LIST.to_vec(){
            let mut rng = SmallRng::seed_from_u64(seed);
            let mut dnfbitmask_generator = RandomDNFBistmaskGenerator::new(&conjunctions, &mut rng);
            let mut formula = dnfbitmask_generator.generate();
            
            let count_before = formula.len();

            let mut dnfbitmask_mutation = SimpleDNFBitmaskMutation::new(&mut rng);
            dnfbitmask_mutation.mutate(&mut formula);

            let count_after = formula.len();

            assert_eq!(count_after.abs_diff(count_before), 1)
        }
    }

    // TO DO: make a test for DNFVecGenerator

    // #[test]
    // fn test_go_terms(){
    //     let seed = 64;
    //     let go: MinimalCsrOntology = get_go(); 
    //     dbg!(go.len());
    //     let mut rng = SmallRng::seed_from_u64(seed);
    //     // let rnd_index = rng.random_range(0..go.len());

    //     // let result_term = go.iter_term_ids().nth(rnd_index);

    //     // if let Some(new_term) = result_term{
    //     //     dbg!(new_term);
    //     //     let children_terms: Vec<&TermId>  = go.iter_child_ids(new_term).collect();
    //     //     dbg!(children_terms.len());
    //     // }

    //     for term in go.iter_term_ids(){
    //         dbg!(go.iter_child_ids(dbg!(term)).count());
    //         dbg!(go.iter_child_ids(term).collect::<Vec<_>>());
    //         dbg!(go.iter_descendant_ids(term).count());
    //         dbg!(go.iter_parent_ids(term).count());
    //         dbg!(go.iter_parent_ids(term).collect::<Vec<_>>());
    //         dbg!(go.iter_ancestor_ids(term).count());
    //         println!("\n\n\n");
    //     }
    // }



}