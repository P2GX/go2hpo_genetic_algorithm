use rand::Rng;
use ontolius::TermId;

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

// Before it was ConjunctionGenerator, but the Disjunctions Generators can implement the same kind of trait, so 
// I renamed it to Formula Generator to be more generic
pub trait FormulaGenerator{
    fn generate(&self) -> Conjunction;
}

pub struct RedundantRandomConjunctionGenerator<'a> {
    n_terms: usize,
    terms: &'a Vec<TermId>,
}

impl<'a> FormulaGenerator for RedundantRandomConjunctionGenerator<'a> {
    fn generate(&self) -> Conjunction {
        let mut terms : Vec<TermObservation> = Vec::new();
        // Select them randomly
        for _ in 0..self.n_terms{
            terms.push(self.randomly_pick());
        }
        Conjunction{ term_observations: terms }
    }
}

impl<'a> RedundantRandomConjunctionGenerator<'a>{
        //randomly generate a TermObeservation
        fn randomly_pick(&self) -> TermObservation{
            //randomly pick a term 
            let mut rng = rand::rng();
            let term_id = self.terms[rng.random_range(0..self.terms.len())].clone();
            //randomly choose whether true or false
            let is_excluded: bool = rng.random(); 
            return TermObservation::new(term_id, is_excluded)
        }
}





#[cfg(test)]
mod tests {
    use super::*;
    use lazy_static::lazy_static;

    lazy_static!{
        static ref small_test_terms:Vec<TermId> = vec!["GO:0051146", "GO:0052693", "GO:0005634"].into_iter().map(|s| s.parse().unwrap()).collect();
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
        assert!(small_test_terms.contains(&observation.term_id), "Term ID is not in the list");
    }


    #[test]
    fn test_generate(){
        let generator = RedundantRandomConjunctionGenerator {
            n_terms: 2,
            terms: &small_test_terms,
        };
        let conjunction: Conjunction = generator.generate();
        println!("{:?}", conjunction);
    }
}
