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



//SELECTION OPERATOR

pub trait Selection<T> {
    fn select(&self, population: &[T]) -> T;
}

