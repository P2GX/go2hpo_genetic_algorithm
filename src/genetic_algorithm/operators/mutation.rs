use crate::genetic_algorithm::Solution;
use rand::prelude::*;
use crate::logical_formula::ConjunctionGenerator;

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

pub struct ConjunctionMutation<O> {
    go: O,
}

impl<O> Mutation<Conjunction> for ConjunctionMutation<O>
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
            7 => self.mutate_with_child_term(formula),
            _ => panic!("A random number outside of the range has been generated. No associated mutation"),
        }
    }
}

impl<O> ConjunctionMutation<O>
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

    /// Toggle the status of a gene expression term from lo to hi or vice versa (And also of GO terms with is_excluded)
    pub fn toggle_term_status(&self, formula: &mut Conjunction) {
        // Is it better to differentiate between GO terms and Gene Expression toggle mutations by
        // creating two separate functions or should it be managed by only one function

        // Toggle go terms
        let mut rng = rand::rng();
        //if this function will manage the toggle of all the different kind of observations, the random number can be
        // obtained from 0 to the length of ALL the annotation terms of the Conjunction
        let rnd_index = rng.random_range(0..formula.term_observations.len());
        let mut term_ob = formula
            .term_observations
            .get_mut(rnd_index)
            .expect("It should return a term");
        term_ob.is_excluded = !term_ob.is_excluded;
        // To do also for the rest
        todo!()
    }

    /// Delete a gene expression term
    pub fn delete_gene_expression_term(&self, formula: &mut Conjunction) {
        todo!()
    }
    /// Add a gene expression term from a random tissue
    pub fn add_gene_expression_term(&self, formula: &mut Conjunction) {
        todo!()
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


pub struct SimpleDNFVecMutation<O, G> 
where 
O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
G: ConjunctionGenerator,
{
    conjunction_mutation: ConjunctionMutation<O>,
    conjunction_generator: G,
}


impl<O, G> Mutation<DNFVec> for SimpleDNFVecMutation<O, G>
where
    O: HierarchyWalks + OntologyTerms<SimpleMinimalTerm>,
    G: ConjunctionGenerator, {
    fn mutate(&self, formula: &mut DNFVec) {
        let mut rng = rand::rng();
        let rnd_num = rng.random_range(0..=7);
        match rnd_num {
            0 => self.mutate_conjunction(formula),
            1 => self.add_random_conjunction(formula),
            2 => self.remove_random_conjunction(formula),
            _ => todo!(),
        }
    }
}


impl<O, G> SimpleDNFVecMutation<O, G>
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
