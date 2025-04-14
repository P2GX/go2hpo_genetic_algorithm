use super::{Conjunction, TermObservation, TissueExpression};
use super::DgeState;
use ontolius::TermId;
use rand::{seq::SliceRandom, Rng};
use std::{collections::{HashMap, HashSet}};

pub trait FormulaGenerator {
    type Output;

    fn generate(&mut self) -> Self::Output;
}


pub trait ConjunctionGenerator: FormulaGenerator<Output = Conjunction> {}


/// Creates Conjunctions randomly that might have redundant terms
pub struct RedundantRandomConjunctionGenerator<'a, R> {
    n_go_terms: usize,
    go_terms: &'a Vec<TermId>,
    n_tissue_terms: usize,
    tissue_terms: &'a Vec<String>,
    rng: R,
}

impl<'a, R> FormulaGenerator  for RedundantRandomConjunctionGenerator<'a, R> 
where
    R: Rng,
{
    type Output = Conjunction;

    fn generate(&mut self) -> Conjunction {
        let mut go_terms: Vec<TermObservation> = Vec::new();
        // Select go terms randomly
        for _ in 0..self.n_go_terms {
            go_terms.push(self.select_random_observation());
        }

        let mut tissue_annots: Vec<TissueExpression> = Vec::new();
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

impl<'a, R> ConjunctionGenerator for RedundantRandomConjunctionGenerator<'a, R> where R: Rng {}

impl<'a, R> RedundantRandomConjunctionGenerator<'a, R> 
where 
    R: Rng,
{
    //randomly generate a TermObeservation
    fn select_random_observation(&mut self) -> TermObservation {
        //randomly pick a term
        let term_id = self.go_terms[self.rng.random_range(0..self.go_terms.len())].clone();
        //randomly choose whether true or false
        let is_excluded: bool = self.rng.random();
        return TermObservation::new(term_id, is_excluded);
    }

    //randomly generate a TissueExpression
    fn select_random_tissue_annot(&mut self) -> TissueExpression {
        let term_id = self.tissue_terms[self.rng.random_range(0..self.tissue_terms.len())].clone();

        let state = DgeState::get_random(&mut self.rng);

        return TissueExpression::new(term_id, state);
    }
}

/// Creates Conjunctions randomly whith non-repeating terms
pub struct RandomConjunctionGenerator<'a, R> {
    n_go_terms: usize,
    go_terms: &'a Vec<TermId>,
    n_tissue_terms: usize,
    tissue_terms: &'a Vec<String>,
    rng: R,
}

impl<'a, R> FormulaGenerator for RandomConjunctionGenerator<'a, R> 
where 
    R: Rng,
{
    type Output = Conjunction;

    fn generate(&mut self) -> Conjunction {

        // Shuffle the annotation vectors
        let mut shuffled_go_terms = self.go_terms.clone();
        shuffled_go_terms.shuffle(&mut self.rng);

        let mut shuffled_tissue_terms = self.tissue_terms.clone();
        shuffled_tissue_terms.shuffle(&mut self.rng);


        // Take the first n according to the value specified in the relative field
        let chosen_go_terms: Vec<TermObservation> = shuffled_go_terms
            .iter()
            .take(self.n_go_terms)
            .map(|term_id| TermObservation::new(term_id.clone(), self.rng.random()))
            .collect();

        let chosen_tissues: Vec<TissueExpression> = shuffled_tissue_terms
            .iter()
            .take(self.n_go_terms)
            .map(|term_id| TissueExpression::new(term_id.clone(), DgeState::get_random(&mut self.rng)))
            .collect();


        Conjunction {
            term_observations: chosen_go_terms,
            tissue_expressions: chosen_tissues,
        }
    }
}

impl<'a, R> ConjunctionGenerator for RandomConjunctionGenerator<'a, R> where R: Rng {}

/// Picks some terms from a list of genes
pub struct GenePickerConjunctionGenerator {
    max_n_go_terms: usize,
    max_n_tissue_annnots: usize,

    //list of genes 
    symbol_to_direct_annotations: HashMap<String, HashSet<TermId>>,
}

impl FormulaGenerator for GenePickerConjunctionGenerator {
    type Output = Conjunction;

    fn generate(&mut self) -> Conjunction {
        todo!();
    }
}

impl ConjunctionGenerator for GenePickerConjunctionGenerator{}

//TO DO: fare dei test per testare anche RandomConjunctionGenerator

#[cfg(test)]
mod tests {
    use super::*;
    use lazy_static::lazy_static;
    use rand::rng;

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
        let mut generator = RedundantRandomConjunctionGenerator {
            n_go_terms: 2,
            go_terms: &small_test_terms,
            n_tissue_terms: 2,
            tissue_terms: &small_test_tissues,
            rng: rng(),
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

    fn _test_generate_random_conjunction<T: ConjunctionGenerator>(conjunction_generator: &mut T){
        let conjunction: Conjunction = conjunction_generator.generate();
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
    fn test_generate_random_conjunction_redundant() {
        let mut generator = RedundantRandomConjunctionGenerator {
            n_go_terms: 2,
            go_terms: &small_test_terms,
            n_tissue_terms: 2,
            tissue_terms: &small_test_tissues,
            rng: rng(), 
        };
        _test_generate_random_conjunction(&mut generator);
    }

    #[test]
    fn test_generate_random_conjunction() {
        let mut generator = RandomConjunctionGenerator {
            n_go_terms: 2,
            go_terms: &small_test_terms,
            n_tissue_terms: 2,
            tissue_terms: &small_test_tissues,
            rng: rng(), 
        };
        _test_generate_random_conjunction(&mut generator);
    }

}
