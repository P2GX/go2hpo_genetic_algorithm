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
    pub fn get_random() -> Self{
        let mut rng = rand::rng();
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
pub struct TermExpression{
    pub term_id: String,
    pub state: DgeState,
}

impl TermExpression {
    pub fn new(term_id: String, state: DgeState) -> Self{
        Self {term_id, state}
    }
}

impl PartialEq for TermExpression {
    fn eq(&self, other: &Self) -> bool {
        // Define custom equality criteria here
        self.term_id == other.term_id && self.state == other.state
    }
}

// impl TermAnnotation<String, DgeState> for TermExpression{
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
    pub tissue_expressions: Vec<TermExpression>,
}

impl PartialEq for Conjunction{
    fn eq(&self, other: &Self) -> bool {
        self.term_observations == other.term_observations && self.tissue_expressions == other.tissue_expressions
    }
}


impl Conjunction{
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




pub trait ConjunctionGenerator {
    fn generate(&self) -> Conjunction;
}

/// Creates Conjunctions randomly that might have redundant terms
pub struct RedundantRandomConjunctionGenerator<'a> {
    n_go_terms: usize,
    go_terms: &'a Vec<TermId>,
    n_tissue_terms: usize,
    tissue_terms: &'a Vec<String>,
}

impl<'a> ConjunctionGenerator for RedundantRandomConjunctionGenerator<'a> {
    fn generate(&self) -> Conjunction {
        let mut go_terms: Vec<TermObservation> = Vec::new();
        // Select go terms randomly
        for _ in 0..self.n_go_terms {
            go_terms.push(self.select_random_observation());
        }

        let mut tissue_annots: Vec<TermExpression> = Vec::new();
        // Select tissues annotations randomly
        for _ in 0..self.n_tissue_terms {
            tissue_annots.push(self.select_random_tissue_annot());
        }
        Conjunction {
            term_observations: go_terms,
            tissue_expressions: tissue_annots,
        }
    }
}

impl<'a> RedundantRandomConjunctionGenerator<'a> {
    //randomly generate a TermObeservation
    fn select_random_observation(&self) -> TermObservation {
        //randomly pick a term
        let mut rng = rand::rng();
        let term_id = self.go_terms[rng.random_range(0..self.go_terms.len())].clone();
        //randomly choose whether true or false
        let is_excluded: bool = rng.random();
        return TermObservation::new(term_id, is_excluded);
    }

    //randomly generate a TermExpression
    fn select_random_tissue_annot(&self) -> TermExpression {
        let mut rng = rand::rng();
        let term_id = self.tissue_terms[rng.random_range(0..self.tissue_terms.len())].clone();

        let state = DgeState::get_random();

        return TermExpression::new(term_id, state);
    }
}

/// Creates Conjunctions randomly whith non-repeating terms
pub struct RandomConjunctionGenerator<'a> {
    n_go_terms: usize,
    go_terms: &'a Vec<TermId>,
    n_tissue_terms: usize,
    tissue_terms: &'a Vec<String>,
}

impl<'a> ConjunctionGenerator for RandomConjunctionGenerator<'a> {
    fn generate(&self) -> Conjunction {
        let mut rng = rand::rng();

        // Shuffle the annotation vectors
        let mut shuffled_go_terms = self.go_terms.clone();
        shuffled_go_terms.shuffle(&mut rng);

        let mut shuffled_tissue_terms = self.tissue_terms.clone();
        shuffled_tissue_terms.shuffle(&mut rng);


        // Take the first n according to the value specified in the relative field
        let chosen_go_terms: Vec<TermObservation> = shuffled_go_terms
            .iter()
            .take(self.n_go_terms)
            .map(|term_id| TermObservation::new(term_id.clone(), rng.random()))
            .collect();

        let chosen_tissues: Vec<TermExpression> = shuffled_tissue_terms
        .iter()
        .take(self.n_go_terms)
        .map(|term_id| TermExpression::new(term_id.clone(), DgeState::get_random()))
        .collect();


        Conjunction {
            term_observations: chosen_go_terms,
            tissue_expressions: chosen_tissues,
        }
    }
}

/// Picks some terms from a list of genes
pub struct GenePickerConjunctionGenerator {
    max_n_go_terms: usize,
    max_n_tissue_annnots: usize,

    //list of genes 
    symbol_to_direct_annotations: HashMap<String, HashSet<TermId>>,
}

impl ConjunctionGenerator for GenePickerConjunctionGenerator {
    fn generate(&self) -> Conjunction {
        todo!();
    }
}

//TO DO: fare dei test per testare anche RandomConjunctionGenerator

#[cfg(test)]
mod tests {
    use super::*;
    use lazy_static::lazy_static;

    lazy_static! {
        static ref small_test_terms: Vec<TermId> = vec!["GO:0051146", "GO:0052693", "GO:0005634"]
            .into_iter()
            .map(|s| s.parse().unwrap())
            .collect();

        static ref small_test_tissues: Vec<String> = vec!["Colon_Transverse_Muscularis".to_string(),
                                                         "Colon_Transverse_Mixed_Cell".to_string(),
                                                          "Colon_Transverse_Muscularis".to_string(),
                                                          "Testis".to_string(),
                                                          "Small_Intestine_Terminal_Ileum_Mixed_Cell".to_string()];
    }

    #[test]
    fn test_randomly_pick() {
        // let small_test_terms:Vec<TermId> = vec!["GO:0051146", "GO:0052693", "GO:0005634"].into_iter().map(|s| s.parse().unwrap()).collect();
        let generator = RedundantRandomConjunctionGenerator {
            n_go_terms: 2,
            go_terms: &small_test_terms,
            n_tissue_terms: 2,
            tissue_terms: &small_test_tissues,
        };

        let observation = generator.select_random_observation();

        // println!("{:?}", observation);
        assert!(
            small_test_terms.contains(&observation.term_id),
            "Term ID is not in the list"
        );

        let tissue = generator.select_random_tissue_annot();
        
        assert!(
            small_test_tissues.contains(&tissue.term_id),
            "Tissue ID is not in the list"
        );

    }

    #[test]
    fn test_generate_redundant_random() {
        let generator = RedundantRandomConjunctionGenerator {
            n_go_terms: 2,
            go_terms: &small_test_terms,
            n_tissue_terms: 2,
            tissue_terms: &small_test_tissues,
        };
        let conjunction: Conjunction = generator.generate();
        println!("{:?}", conjunction);

        assert_eq!(conjunction.len(), 4, "The real number of terms differs from the expected one");
        
        for term_obs in conjunction.term_observations.iter(){
            assert!(
                small_test_terms.contains(&term_obs.term_id),
                "Term ID is not in the list"
            );
        }

        for tissue_term in conjunction.tissue_expressions{
            assert!(
                small_test_tissues.contains(&tissue_term.term_id),
                "Term ID is not in the list"
            );
        }

    }

    #[test]
    fn test_generate_random() {
        let generator = RedundantRandomConjunctionGenerator {
            n_go_terms: 2,
            go_terms: &small_test_terms,
            n_tissue_terms: 2,
            tissue_terms: &small_test_tissues,
        };
        let conjunction: Conjunction = generator.generate();
        println!("{:?}", conjunction);

        assert_eq!(conjunction.len(), 4, "The real number of terms differs from the expected one");
        
        for term_obs in conjunction.term_observations{
            assert!(
                small_test_terms.contains(&term_obs.term_id),
                "Term ID is not in the list"
            );
        }

        for tissue_term in conjunction.tissue_expressions{
            assert!(
                small_test_tissues.contains(&tissue_term.term_id),
                "Term ID is not in the list"
            );
        }

    }

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
            } else if let Some(v) = items.downcast_ref::<Vec<TermExpression>>() {
                println!("Tissue Expressions count = {:?}", v.len());
            } else {
                println!("Unknown Type");
            }
        }


    }
}
