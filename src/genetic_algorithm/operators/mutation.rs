use crate::genetic_algorithm::Solution;
use rand::prelude::*;
use crate::logical_formula::{ConjunctionGenerator, DgeState, TermExpression};
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
    fn mutate(&self, formula: &mut T);
}

pub struct ConjunctionMutation<'a, O> {
    go: &'a O,
    gtex: &'a GtexSummary,
}

impl<O> Mutation<Conjunction> for ConjunctionMutation<'_, O>
where
    O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
{
    fn mutate(&self, formula: &mut Conjunction) {
        let mut rng = rand::rng();
        let rnd_num = rng.random_range(0..=7);
        match rnd_num {
            0 => self.mutate_with_parent_term(formula),
            1 => self.mutate_with_child_term(formula),
            2 => self.delete_random_term(formula),
            3 => self.add_random_term(formula),
            4 => self.toggle_term_status(formula),
            5 => self.delete_gene_expression_term(formula),
            6 => self.add_gene_expression_term(formula),
            7 => self.toggle_gene_expression_term(formula),
            _ => panic!("A random number outside of the range has been generated. No associated mutation"),
        }
    }
}

impl<O> ConjunctionMutation<'_, O>
where
    O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
{
    /// Exchange an HPO term with a parent
    pub fn mutate_with_parent_term(&self, formula: &mut Conjunction) {
        // Get a term from a random index. Maybe a more sophisticated way will be used in the future
        let mut rng = rand::rng();
        let rnd_index = rng.random_range(0..formula.term_observations.len());
        let mut term_ob = formula
            .term_observations
            .get_mut(rnd_index)
            .expect("It should return a term");

        // Get the parents' term IDs, I collect them to know the length without consuming it
        let parents_ids: Vec<&TermId> = self.go.iter_parent_ids(&term_ob.term_id).collect();

        // Select one of them randomly and subsitute it with the current termId in the formula
        let rnd_index = rng.random_range(0..parents_ids.len());
        term_ob.term_id = parents_ids.get(rnd_index).copied().unwrap().clone();
    }

    ///  Exchange an HPO term with one of its children
    pub fn mutate_with_child_term(&self, formula: &mut Conjunction) {
        // Get a term from a random index. Maybe a more sophisticated way will be used in the future
        let mut rng = rand::rng();
        let rnd_index = rng.random_range(0..formula.term_observations.len());
        let mut term_ob = formula
            .term_observations
            .get_mut(rnd_index)
            .expect("It should return a term");

        // Get the parents' term IDs, I collect them to know the length without consuming it
        let child_ids: Vec<&TermId> = self.go.iter_child_ids(&term_ob.term_id).collect();

        // Select one of them randomly and subsitute it with the current termId in the formula
        let rnd_index = rng.random_range(0..child_ids.len());
        term_ob.term_id = child_ids.get(rnd_index).copied().unwrap().clone();
    }

    /// Delete an HPO term
    pub fn delete_random_term(&self, formula: &mut Conjunction) {
        // Get a term from a random index. Maybe a more sophisticated way will be used in the future
        let mut rng = rand::rng();
        let rnd_index = rng.random_range(0..formula.term_observations.len());
        formula.term_observations.remove(rnd_index);
    }

    /// Add a random HPO term
    pub fn add_random_term(&self, formula: &mut Conjunction) {
        let mut rng = rand::rng();
        let rnd_index = rng.random_range(0..self.go.len());
        if let Some(new_term) = self.go.iter_term_ids().nth(rnd_index) {
            let term_obs = TermObservation::new(new_term.clone(), rng.random_bool(0.5));
            formula.term_observations.push(term_obs);
        }
    }

    /// Toggle the status of GO terms with is_excluded
    pub fn toggle_term_status(&self, formula: &mut Conjunction) {

        let mut rng = rand::rng();

        let rnd_index = rng.random_range(0..formula.term_observations.len());
        let mut term_ob = formula
            .term_observations
            .get_mut(rnd_index)
            .expect("It should return a term");
        term_ob.is_excluded = !term_ob.is_excluded;
    }

    /// Delete a gene expression term
    pub fn delete_gene_expression_term(&self, formula: &mut Conjunction) {
        let mut rng = rand::rng();
        let rnd_index = rng.random_range(0..formula.tissue_expressions.len());
        formula.tissue_expressions.remove(rnd_index);
    }
    
    /// Add a gene expression term from a random tissue
    pub fn add_gene_expression_term(&self, formula: &mut Conjunction) {
        let mut rng = rand::rng();
        let tissues = self.gtex.metadata.get_tissue_names();
        
        let rnd_index = rng.random_range(0..tissues.len());
        
        if let Some(tissue_name) = tissues.get(rnd_index){
            let tissue_expr = TermExpression::new(tissue_name.clone(), DgeState::get_random());
            formula.tissue_expressions.push(tissue_expr)
        }
    }

    // toggle of a gene expression term from lo to hi or vice versa
    pub fn toggle_gene_expression_term(&self, formula: &mut Conjunction) {
        let mut rng = rand::rng();
        let rnd_index = rng.random_range(0..formula.tissue_expressions.len());
        let mut tissue_term = formula.tissue_expressions.get_mut(rnd_index).expect("It should return a TermExpression");
        
        let dge_states = [DgeState::Down, DgeState::Normal, DgeState::Up];

        let new_state = match tissue_term.state {
            DgeState::Down =>  [DgeState::Normal, DgeState::Up].choose(&mut rng).expect("Should return one state"),
            DgeState::Normal =>  [DgeState::Down, DgeState::Up].choose(&mut rng).expect("Should return one state"),
            DgeState::Up =>  [DgeState::Normal, DgeState::Down].choose(&mut rng).expect("Should return one state"),
        };

        tissue_term.state =  new_state.clone();
        
    }

    
}

pub struct SimpleDNFBitmaskMutation;

impl Mutation<DNFBitmask<'_>> for SimpleDNFBitmaskMutation {
    fn mutate(&self, formula: &mut DNFBitmask) {
        let n_conjunctions = formula.total_conjunctions_count();
        let mut rng = rand::rng();
        let rnd_index = rng.random_range(0..n_conjunctions);

        formula
            .toggle_conjunction(rnd_index)
            .expect("It should toggle the conjunction");
    }
}


pub struct SimpleDNFVecMutation<'a, O, G> 
where 
O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
G: ConjunctionGenerator,
{
    conjunction_mutation: ConjunctionMutation<'a, O>,
    conjunction_generator: G,
}


impl<'a, O, G> Mutation<DNFVec> for SimpleDNFVecMutation<'a, O, G>
where
    O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
    G: ConjunctionGenerator, {
    fn mutate(&self, formula: &mut DNFVec) {
        let mut rng = rand::rng();
        let rnd_num = rng.random_range(0..=2);
        match rnd_num {
            0 => self.mutate_conjunction(formula),
            1 => self.add_random_conjunction(formula),
            2 => self.remove_random_conjunction(formula),
            _ => panic!("A random number outside of the range has been generated. No associated mutation"),
        }
    }
}


impl<'a, O, G> SimpleDNFVecMutation<'a, O, G>
where
    O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
    G: ConjunctionGenerator,
{
    /// Choose one conjunction randomly and then modifies it with ConjunctionMutation
    pub fn mutate_conjunction(&self, formula: &mut DNFVec) {
        let mut conjunctions = formula.get_mut_active_conjunctions();
        let mut rng = rand::rng();
        let rnd_index = rng.random_range(0..conjunctions.len());
        let mut conjunction = conjunctions
            .get_mut(rnd_index)
            .expect("It should return a conjunction");
        self.conjunction_mutation.mutate(&mut conjunction);
    }

    /// creates a new conjunction and assigns it to the formula
    pub fn add_random_conjunction(&self, formula: &mut DNFVec) {
        // In a future implementatio the conjunction might be copied from one of the already existing and best performing
        let conjunction: Conjunction = self.conjunction_generator.generate();
        formula.activate_conjunction(conjunction).expect("The conjunction should be added without errors")
    }

    /// Removes a conjunction randommly
    pub fn remove_random_conjunction(&self, formula: &mut DNFVec){
        //In a next implementation the probability of being removed might depend on the conjunction performance
        let mut conjunctions = formula.get_mut_active_conjunctions();
        let mut rng = rand::rng();
        let rnd_index = rng.random_range(0..conjunctions.len());
        conjunctions.remove(rnd_index);
    }
}

pub struct BiasedDNFMutation;

impl Mutation<DNFBitmask<'_>> for BiasedDNFMutation {
    fn mutate(&self, formula: &mut DNFBitmask) {
        todo!()
    }
}
