use ontolius::TermId;
use rand::{seq::SliceRandom, Rng};

#[derive(Debug)]
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

#[derive(Debug)]
pub struct Conjunction {
    pub term_observations: Vec<TermObservation>,
}

impl Conjunction{
    pub fn len(&self) -> usize{
        return self.term_observations.len() ;
    }
}

// Before it was ConjunctionGenerator, but the Disjunctions Generators can implement the same kind of trait, so
// I renamed it to Formula Generator to be more generic
pub trait FormulaGenerator {
    fn generate(&self) -> Conjunction;
}

//Creates Conjunctions randomly that might have redundant terms
pub struct RedundantRandomConjunctionGenerator<'a> {
    n_terms: usize,
    terms: &'a Vec<TermId>,
}

impl<'a> FormulaGenerator for RedundantRandomConjunctionGenerator<'a> {
    fn generate(&self) -> Conjunction {
        let mut terms: Vec<TermObservation> = Vec::new();
        // Select them randomly
        for _ in 0..self.n_terms {
            terms.push(self.randomly_pick());
        }
        Conjunction {
            term_observations: terms,
        }
    }
}

impl<'a> RedundantRandomConjunctionGenerator<'a> {
    //randomly generate a TermObeservation
    fn randomly_pick(&self) -> TermObservation {
        //randomly pick a term
        let mut rng = rand::rng();
        let term_id = self.terms[rng.random_range(0..self.terms.len())].clone();
        //randomly choose whether true or false
        let is_excluded: bool = rng.random();
        return TermObservation::new(term_id, is_excluded);
    }
}

//Creates Conjunctions randomly whith non-repeating terms
pub struct RandomConjunctionGenerator<'a> {
    n_terms: usize,
    terms: &'a Vec<TermId>,
}

impl<'a> FormulaGenerator for RandomConjunctionGenerator<'a> {
    fn generate(&self) -> Conjunction {
        let mut rng = rand::rng();

        let mut shuffled_terms = self.terms.clone();
        shuffled_terms.shuffle(&mut rng);

        let selected_terms: Vec<TermObservation> = shuffled_terms
            .iter()
            .take(self.n_terms)
            .map(|term_id| TermObservation::new(term_id.clone(), rng.random()))
            .collect();

        Conjunction {
            term_observations: selected_terms,
        }
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
    }

    #[test]
    fn test_randomly_pick() {
        // let small_test_terms:Vec<TermId> = vec!["GO:0051146", "GO:0052693", "GO:0005634"].into_iter().map(|s| s.parse().unwrap()).collect();
        let generator = RedundantRandomConjunctionGenerator {
            n_terms: 2,
            terms: &small_test_terms,
        };

        let observation = generator.randomly_pick();

        // println!("{:?}", observation);
        assert!(
            small_test_terms.contains(&observation.term_id),
            "Term ID is not in the list"
        );
    }

    #[test]
    fn test_generate_redundant_random() {
        let n_terms = 2;
        let generator = RedundantRandomConjunctionGenerator {
            n_terms: 2,
            terms: &small_test_terms,
        };
        let conjunction: Conjunction = generator.generate();
        println!("{:?}", conjunction);

        assert_eq!(conjunction.len(), n_terms, "The real number of terms differs from the expected one");
        
        for term_obs in conjunction.term_observations{
            assert!(
                small_test_terms.contains(&term_obs.term_id),
                "Term ID is not in the list"
            );
        }

    }

    #[test]
    fn test_generate_random() {
        let n_terms = 2;
        let generator = RedundantRandomConjunctionGenerator {
            n_terms: 2,
            terms: &small_test_terms,
        };
        let conjunction: Conjunction = generator.generate();
        println!("{:?}", conjunction);

        assert_eq!(conjunction.len(), n_terms, "The real number of terms differs from the expected one");
        
        for term_obs in conjunction.term_observations{
            assert!(
                small_test_terms.contains(&term_obs.term_id),
                "Term ID is not in the list"
            );
        }

    }

    // #[test]
    // fn test_generate(){
        
    // }
}
