//! Core literal types (GO term observations, tissue expressions) and Conjunction.
use ontolius::TermId;
use rand::{seq::SliceRandom, Rng};
use serde::{Deserialize, Serialize};
use std::fmt;
use std::{
    any::Any,
    collections::{HashMap, HashSet},
};

use crate::logical_formula::base::Formula;

pub trait TermAnnotation {}
// T should be TermId for TermObservation (GO) and a String (at the moment, maybe something more complicated in the future) for Gtex Expression Data
// V should be a bool for TermObservation (GO) and an Enum {UP, DOWN, NORMAL} for GtexExpressionData
// pub trait TermAnnotation<T, V>{
//     fn get_term(&self) -> &T;
//     fn is_annotated_as(&self, annotation_type: V) -> bool;
// }

/// GO term literal used inside a conjunction. `is_excluded=true` means NOT(term).
#[derive(Debug, Clone)]
pub struct TermObservation {
    /// GO term identifier.
    pub term_id: TermId,
    /// Whether the term is negated in the conjunction.
    pub is_excluded: bool,
}

impl TermObservation {
    /// Create a GO literal; `is_excluded=true` means NOT(term).
    /// Inclusion matches direct term or its descendants (checked elsewhere).
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

impl fmt::Display for TermObservation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_excluded {
            write!(f, "NOT({})", self.term_id)
        } else {
            write!(f, "{}", self.term_id)
        }
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

/// Differential gene-expression state (up/down/normal).
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DgeState {
    Up,
    Down,
    Normal,
}

impl DgeState {
    /// Randomly draw any of Up/Down/Normal.
    pub fn get_random<R: Rng>(rng: &mut R) -> Self {
        // let mut rng = rand::rng();
        let rnd_state_i = rng.random_range(0..3);
        match rnd_state_i {
            0 => DgeState::Down,
            1 => DgeState::Normal,
            2 => DgeState::Up,
            _ => panic!("Number out of the range has been generated for DgeState"),
        }
    }

    /// Returns only UP or DOWN, excluding NORMAL.
    pub fn get_random_up_down_only<R: Rng>(rng: &mut R) -> Self {
        let rnd_state_i = rng.random_range(0..2);
        match rnd_state_i {
            0 => DgeState::Down,
            1 => DgeState::Up,
            _ => panic!(
                "Number out of the range has been generated for DgeState::get_random_up_down_only"
            ),
        }
    }
}

pub type TissueId = String;

/// Tissue expression literal paired with a DGE state.
#[derive(Debug, Clone, Eq, Hash, Serialize, Deserialize)]
pub struct TissueExpression {
    pub term_id: TissueId,
    pub state: DgeState,
}

impl TissueExpression {
    /// Create a tissue expression literal.
    pub fn new(term_id: TissueId, state: DgeState) -> Self {
        Self { term_id, state }
    }

    pub fn into_down(&self) -> Self {
        TissueExpression {
            term_id: self.term_id.clone(),
            state: DgeState::Down,
        }
    }

    pub fn into_up(&self) -> Self {
        TissueExpression {
            term_id: self.term_id.clone(),
            state: DgeState::Up,
        }
    }
}

impl PartialEq for TissueExpression {
    fn eq(&self, other: &Self) -> bool {
        // Define custom equality criteria here
        self.term_id == other.term_id && self.state == other.state
    }
}

impl fmt::Display for DgeState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DgeState::Up => write!(f, "UP"),
            DgeState::Down => write!(f, "DOWN"),
            DgeState::Normal => write!(f, "NORMAL"),
        }
    }
}

impl fmt::Display for TissueExpression {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.state {
            DgeState::Up => write!(f, "UP({})", self.term_id),
            DgeState::Down => write!(f, "DOWN({})", self.term_id),
            DgeState::Normal => write!(f, "NORMAL({})", self.term_id),
        }
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

/// A conjunction (logical AND) of GO term observations and tissue expressions.
/// Example display: `(GO:0008150 AND NOT(GO:0003674) AND UP(Liver))`.
#[derive(Debug, Clone)]
pub struct Conjunction {
    /// GO term literals (included or excluded).
    pub term_observations: Vec<TermObservation>,
    /// Tissue expression literals (UP/DOWN/NORMAL).
    pub tissue_expressions: Vec<TissueExpression>,
}

impl PartialEq for Conjunction {
    fn eq(&self, other: &Self) -> bool {
        self.term_observations == other.term_observations
            && self.tissue_expressions == other.tissue_expressions
    }
}

impl Conjunction {
    /// Empty conjunction.
    pub fn new() -> Self {
        Self {
            term_observations: vec![],
            tissue_expressions: vec![],
        }
    }

    /// Build from provided literal vectors.
    pub fn from(
        term_observations: Vec<TermObservation>,
        tissue_expressions: Vec<TissueExpression>,
    ) -> Self {
        Self {
            term_observations,
            tissue_expressions,
        }
    }

    /// Iterate over fields without naming (term_observations then tissue_expressions).
    pub fn iter(&self) -> ConjunctionIterator {
        ConjunctionIterator::new(self)
    }

    /// Iterate over fields with their names (used by generic traversal).
    pub fn named_iter(&self) -> ConjunctionNamedIterator {
        ConjunctionNamedIterator::new(self)
    }
}

impl Formula for Conjunction {
    // Sum of all the annotation vectors
    /// Total length of all the annotation terms present in the conjunction
    fn len(&self) -> usize {
        return self.term_observations.len() + self.tissue_expressions.len();
    }
}

impl fmt::Display for Conjunction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Collect all parts (TermObservations + TissueExpressions)
        let mut parts: Vec<String> = Vec::new();

        parts.extend(self.term_observations.iter().map(|t| t.to_string()));
        parts.extend(self.tissue_expressions.iter().map(|t| t.to_string()));

        write!(f, "({})", parts.join(" AND "))
    }
}

pub struct ConjunctionNamedIterator<'a> {
    state: usize,
    conjunction: &'a Conjunction,
}

impl<'a> ConjunctionNamedIterator<'a> {
    /// Iterator that yields ("term_observations", Vec<TermObservation>) then ("tissue_expressions", Vec<TissueExpression>).
    pub fn new(conjunction: &'a Conjunction) -> Self {
        Self {
            state: 0,
            conjunction: conjunction,
        }
    }
}

impl<'a> Iterator for ConjunctionNamedIterator<'a> {
    type Item = (&'a str, &'a dyn Any);

    fn next(&mut self) -> Option<Self::Item> {
        let result = match self.state {
            0 => Some((
                "term_observations",
                &self.conjunction.term_observations as &dyn Any,
            )),
            1 => Some((
                "tissue_expressions",
                &self.conjunction.tissue_expressions as &dyn Any,
            )),
            _ => None,
        };

        self.state += 1;
        result
    }
}

pub struct ConjunctionIterator<'a> {
    state: usize,
    conjunction: &'a Conjunction,
}

impl<'a> ConjunctionIterator<'a> {
    /// Iterator that yields the two fields as `&dyn Any`, in order.
    pub fn new(conjunction: &'a Conjunction) -> Self {
        Self {
            state: 0,
            conjunction: conjunction,
        }
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

    fn test_generate() {
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

    #[test]
    fn test_display_termobservation() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let obs1 = TermObservation::new(t1.clone(), false);
        let obs2 = TermObservation::new(t1, true);

        // println!("{}", obs1); // prints: GO:0051146
        // println!("{}", obs2); // prints: NOT(GO:0051146)

        assert_eq!(format!("{}", obs1), "GO:0051146");
        assert_eq!(format!("{}", obs2), "NOT(GO:0051146)");
    }

    #[test]
    fn test_display_tissueexpression() {
        let t1 = TissueExpression::new("Liver".to_string(), DgeState::Up);
        let t2 = TissueExpression::new("Heart".to_string(), DgeState::Down);
        let t3 = TissueExpression::new("Brain".to_string(), DgeState::Normal);

        assert_eq!(format!("{}", t1), "UP(Liver)");
        assert_eq!(format!("{}", t2), "DOWN(Heart)");
        assert_eq!(format!("{}", t3), "NORMAL(Brain)");
    }

    #[test]
    fn test_display_conjunction() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0051216".parse().unwrap();

        let conj = Conjunction {
            term_observations: vec![
                TermObservation::new(t1, false),
                TermObservation::new(t2, true),
            ],
            tissue_expressions: vec![TissueExpression::new("Liver".to_string(), DgeState::Up)],
        };
        println!("{}", conj);
        assert_eq!(
            format!("{}", conj),
            "(GO:0051146 AND NOT(GO:0051216) AND UP(Liver))"
        );
    }
}
