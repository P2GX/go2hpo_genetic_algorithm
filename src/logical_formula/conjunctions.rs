use ontolius::TermId;
use rand::{seq::SliceRandom, Rng};
use std::{any::Any, collections::{HashMap, HashSet}};





pub trait TermAnnotation{}
// T should be TermId for TermObservation (GO) and a String (at the moment, maybe something more complicated in the future) for Gtex Expression Data
// V should be a bool for TermObservation (GO) and an Enum {UP, DOWN, NORMAL} for GtexExpressionData
// pub trait TermAnnotation<T, V>{
//     fn get_term(&self) -> &T;
//     fn is_annotated_as(&self, annotation_type: V) -> bool; 
// }

#[derive(Debug, Clone)]
pub struct TermObservation {
    pub term_id: TermId,
    pub is_excluded: bool,
}

impl TermObservation {
    pub fn new(term_id: TermId, is_excluded: bool) -> Self {
        Self {
            term_id,
            is_excluded,
        }
    }
}

impl PartialEq for TermObservation {
    fn eq(&self, other: &Self) -> bool {
        // Define custom equality criteria here
        self.term_id == other.term_id && self.is_excluded == other.is_excluded
    }
}


// impl TermAnnotation<TermId, bool> for TermObservation{
//     fn get_term(&self) -> &TermId {
//         return &self.term_id;
//     }
//     fn is_annotated_as(&self, is_excluded: bool) -> bool {
//         return self.is_excluded == is_excluded;
//     }
// }

#[derive(Debug, Clone, PartialEq)]
pub enum DgeState{
    Up,
    Down,
    Normal
}

impl DgeState{
    pub fn get_random<R: Rng>(rng: &mut R) -> Self{
        // let mut rng = rand::rng();
        let rnd_state_i  = rng.random_range(0..3);
        match  rnd_state_i{
            0 => DgeState::Down,
            1 => DgeState::Normal,
            2 => DgeState::Up,
            _ => panic!("Number out of the range has been generated for DgeState"),
        }
    }
}

#[derive(Debug, Clone)]
pub struct TissueExpression{
    pub term_id: String,
    pub state: DgeState,
}

impl TissueExpression {
    pub fn new(term_id: String, state: DgeState) -> Self{
        Self {term_id, state}
    }
}

impl PartialEq for TissueExpression {
    fn eq(&self, other: &Self) -> bool {
        // Define custom equality criteria here
        self.term_id == other.term_id && self.state == other.state
    }
}

// impl TermAnnotation<String, DgeState> for TissueExpression{
//     fn get_term(&self) -> &String {
//         return &self.term_id;
//     }
//     fn is_annotated_as(&self, annotation_type: DgeState) -> bool {
//         return self.state == annotation_type;
//     }
// }


#[derive(Debug, Clone)]
pub struct Conjunction {
    pub term_observations: Vec<TermObservation>,
    pub tissue_expressions: Vec<TissueExpression>,
}

impl PartialEq for Conjunction{
    fn eq(&self, other: &Self) -> bool {
        self.term_observations == other.term_observations && self.tissue_expressions == other.tissue_expressions
    }
}


impl Conjunction{
    
    pub fn new() -> Self{
        Self { term_observations: vec![], tissue_expressions: vec![] }
    }
    
    // Sum of all the annotation vectors
    /// Total length of all the annotation terms present in the conjunction
    pub fn len(&self) -> usize{
        return self.term_observations.len() + self.tissue_expressions.len();
    }

    pub fn iter(&self) -> ConjunctionIterator{
        ConjunctionIterator::new(self)
    }

    pub fn named_iter(&self) -> ConjunctionNamedIterator{
        ConjunctionNamedIterator::new(self)
    }
}

pub struct ConjunctionNamedIterator<'a>{
    state: usize,
    conjunction: &'a Conjunction,
}

impl<'a> ConjunctionNamedIterator<'a>{
    pub fn new(conjunction: &'a Conjunction) -> Self{
        Self{state: 0, conjunction: conjunction}
    }
}

impl<'a> Iterator for ConjunctionNamedIterator<'a> {
    type Item = (&'a str, &'a dyn Any);

    fn next(&mut self) -> Option<Self::Item> {
        let result = match self.state {
            0 => Some(("term_observations", &self.conjunction.term_observations as &dyn Any)),
            1 => Some(("tissue_expressions", &self.conjunction.tissue_expressions as &dyn Any)),
            _ => None,
        };

        self.state += 1;
        result
    }
}


pub struct ConjunctionIterator<'a>{
    state: usize,
    conjunction: &'a Conjunction,
}

impl<'a> ConjunctionIterator<'a>{
    pub fn new(conjunction: &'a Conjunction) -> Self{
        Self{state: 0, conjunction: conjunction}
    }
}

impl<'a> Iterator for ConjunctionIterator<'a> {
    type Item = &'a dyn Any;

    fn next(&mut self) -> Option<Self::Item> {
        let result = match self.state {
            0 => Some(&self.conjunction.term_observations as &dyn Any),
            1 => Some(&self.conjunction.tissue_expressions as &dyn Any),
            _ => None,
        };

        self.state += 1;
        result
    }
}



#[cfg(test)]
mod tests {
    use super::*;
    use lazy_static::lazy_static;
    use rand::rng;

    #[test]
    fn test_generate(){
        let conj = Conjunction {
            term_observations: vec![],
            tissue_expressions: vec![],
        };
    
        let mut iter = ConjunctionNamedIterator::new(&conj);
    
        while let Some((field_name, items)) = iter.next() {
            print!("{}: ", field_name);
            if let Some(v) = items.downcast_ref::<Vec<TermObservation>>() {
                println!("Term Observations count = {:?}", v.len());
            } else if let Some(v) = items.downcast_ref::<Vec<TissueExpression>>() {
                println!("Tissue Expressions count = {:?}", v.len());
            } else {
                println!("Unknown Type");
            }
        }
    }
}
