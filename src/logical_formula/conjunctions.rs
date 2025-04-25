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

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
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

pub type TissueId = String;

#[derive(Debug, Clone, Eq, Hash)]
pub struct TissueExpression{
    pub term_id: TissueId,
    pub state: DgeState,
}

impl TissueExpression {
    pub fn new(term_id: TissueId, state: DgeState) -> Self{
        Self {term_id, state}
    }

    pub fn into_down(&self) -> Self{
        TissueExpression { term_id: self.term_id.clone(), state: DgeState::Down }
    }

    pub fn into_up(&self) -> Self{
        TissueExpression { term_id: self.term_id.clone(), state: DgeState::Up }
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

    #[test]
    fn test_new() {
        let expr = TissueExpression::new("TISSUE1".to_string(), DgeState::Up);
        assert_eq!(expr.term_id, "TISSUE1");
        assert_eq!(expr.state, DgeState::Up);
    }

    #[test]
    fn test_into_down() {
        let expr = TissueExpression::new("TISSUE1".to_string(), DgeState::Up);
        let down_expr = expr.into_down();
        assert_eq!(down_expr.term_id, "TISSUE1");
        assert_eq!(down_expr.state, DgeState::Down);
    }

    #[test]
    fn test_into_up() {
        let expr = TissueExpression::new("TISSUE1".to_string(), DgeState::Down);
        let up_expr = expr.into_up();
        assert_eq!(up_expr.term_id, "TISSUE1");
        assert_eq!(up_expr.state, DgeState::Up);
    }

    #[test]
    fn test_equality() {
        let expr1 = TissueExpression::new("TISSUE1".to_string(), DgeState::Up);
        let expr2 = TissueExpression::new("TISSUE1".to_string(), DgeState::Up);
        let expr3 = TissueExpression::new("TISSUE1".to_string(), DgeState::Down);
        assert_eq!(expr1, expr2);
        assert_ne!(expr1, expr3);
    }

    #[test]
    fn test_new_term_observation() {
        let term_id: TermId = "GO:0051146".parse().unwrap();
        let obs = TermObservation::new(term_id.clone(), true);
        assert_eq!(obs.term_id, term_id);
        assert!(obs.is_excluded);
    }

    #[test]
    fn test_term_observation_equality() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0051146".parse().unwrap();
        let t3: TermId = "GO:0008150".parse().unwrap();

        let obs1 = TermObservation::new(t1, false);
        let obs2 = TermObservation::new(t2, false);
        let obs3 = TermObservation::new(t3, true);

        assert_eq!(obs1, obs2);
        assert_ne!(obs1, obs3);
    }

    #[test]
    fn test_new_conjunction_is_empty() {
        let conj = Conjunction::new();
        assert!(conj.term_observations.is_empty());
        assert!(conj.tissue_expressions.is_empty());
        assert_eq!(conj.len(), 0);
    }

    #[test]
    fn test_conjunction_equality() {
        let t1: TermId = "GO:0008150".parse().unwrap();
        let t2: TermId = "GO:0051146".parse().unwrap();

        let term_obs1 = TermObservation::new(t1.clone(), false);
        let term_obs2 = TermObservation::new(t2.clone(), true);

        let expr1 = TissueExpression::new("TISSUE1".to_string(), DgeState::Up);
        let expr2 = TissueExpression::new("TISSUE2".to_string(), DgeState::Down);
        let expr3 = TissueExpression::new("TISSUE3".to_string(), DgeState::Down);

        let mut conj1 = Conjunction::new();
        conj1.term_observations.push(term_obs1.clone());
        conj1.term_observations.push(term_obs2.clone());
        conj1.tissue_expressions.push(expr1.clone());
        conj1.tissue_expressions.push(expr2.clone());

        let mut conj2 = Conjunction::new();
        conj2.term_observations.push(term_obs1.clone());
        conj2.term_observations.push(term_obs2.clone());
        conj2.tissue_expressions.push(expr1.clone());
        conj2.tissue_expressions.push(expr2.clone());

        let mut conj3 = Conjunction::new();
        conj3.term_observations.push(term_obs1.clone());
        conj3.term_observations.push(term_obs2.clone());
        conj3.tissue_expressions.push(expr1.clone());
        conj3.tissue_expressions.push(expr3.clone());

        assert_eq!(conj1, conj2);
        assert_eq!(conj1.len(), 4);
        assert_ne!(conj1, conj3);
    }
}
