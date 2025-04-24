use lazy_static::lazy_static;

use std::{
    cell::OnceCell,
    collections::{HashMap, HashSet},
};

use ontolius::{ontology::HierarchyWalks, TermId};

use crate::annotations::{GeneAnnotations, GeneSetAnnotations};

use super::{Conjunction, TissueExpression};

pub trait SatisfactionChecker {
    fn is_satisfied(&self, symbol: &str, conjunction: &Conjunction) -> bool;
    fn all_satisfactions(&self, conjunction: &Conjunction) -> HashMap<String, bool>;
}


pub struct NaiveSatisfactionChecker<O> {
    go: O,
    // symbol_to_direct_annotations: HashMap<String, HashSet<TermId>>,
    gene_set_annotations: &'static GeneSetAnnotations, 
}
impl<O> SatisfactionChecker for NaiveSatisfactionChecker<O>
where
    O: HierarchyWalks,
{
    /// Checks if a gene satisfies the annotations of a conjunction
    fn is_satisfied(&self, symbol: &str, conjunction: &Conjunction) -> bool {
        let gene_annotations_lookup = self.gene_set_annotations.get_gene_annotations(symbol);
        match gene_annotations_lookup{
            Some(gene_annotations) => {
                if self.are_term_annotations_satisfied(gene_annotations, conjunction) == false {return false;}
                if self.are_tissue_expressions_satisfied(gene_annotations, conjunction) == false {return false;}
                return true;
            },
            None => panic!("Couldn't find the gene ID in the Gene Annotation Set"),
        }
    }

    /// Checks satisfaction for all genes wrt to a conjunction and returns a HashMap of gene symbols 
    /// and boolean value indicating the satisfaction
    fn all_satisfactions(&self, conjunction: &Conjunction) -> HashMap<String, bool> {
        self.gene_set_annotations
            .get_gene_annotations_map()
            .iter()
            .map(|(symbol, gene_annotations)| {
                let satisfied = self.are_term_annotations_satisfied(gene_annotations, conjunction)
                    && self.are_tissue_expressions_satisfied(gene_annotations, conjunction);
                (symbol.clone(), satisfied)
            })
            .collect()
    }
}

impl<O> NaiveSatisfactionChecker<O> 
where
    O: HierarchyWalks,{
    pub fn new(go: O, gene_set_annotations: &'static GeneSetAnnotations) -> Self {
        Self {
            go,
            gene_set_annotations,
        }
    }


    pub fn are_term_annotations_satisfied(&self, gene_annotations: &GeneAnnotations, conjunction: &Conjunction) -> bool{
        let direct_annotations = gene_annotations.get_term_annotations();
        for term_ob in &conjunction.term_observations {
            //LONGER VERSION:
            match term_ob.is_excluded {
                true => {
                    if direct_annotations.contains(&term_ob.term_id)
                        || self
                            .go
                            .iter_descendant_ids(&term_ob.term_id)
                            .any(|desc| direct_annotations.contains(desc))
                    {
                        // Not satisfied because the gene has a direct or undirect annotation to term_ob, which is excluded.
                        return false;
                    }
                }
                false => {
                    // self.go.
                    if !direct_annotations.contains(&term_ob.term_id)
                        && !self
                            .go
                            .iter_descendant_ids(&term_ob.term_id)
                            .any(|desc| direct_annotations.contains(desc))
                    {
                        return false;
                    }
                }
            }
            // I can change the match with a shorter version (see notes)
        }
        true
    }

    pub fn are_tissue_expressions_satisfied(&self, gene_annotations: &GeneAnnotations, conjunction: &Conjunction) -> bool{
        conjunction.tissue_expressions
            .iter()
            .all(|tissue_expr| gene_annotations.contains_tissue_expressions(tissue_expr))
    }

}

#[cfg(test)]
mod tests {
    use std::{fs::File, hash::Hash, io::BufReader};

    use flate2::bufread::GzDecoder;
    use ontolius::{io::OntologyLoaderBuilder, ontology::csr::MinimalCsrOntology};

    use crate::logical_formula::TermObservation;

    use super::*;

    lazy_static! {
        static ref checker : NaiveSatisfactionChecker<MinimalCsrOntology> = initialize_data(); // NaiveSatisfactionChecker::new(go, map);
    }

    fn initialize_data() -> NaiveSatisfactionChecker<MinimalCsrOntology> {
        let go_path = "data/go/go.toy.json.gz";
        let reader = GzDecoder::new(BufReader::new(
            File::open(go_path).expect("The file should be in the repo"),
        ));

        let loader = OntologyLoaderBuilder::new().obographs_parser().build();
        let go: MinimalCsrOntology = loader
            .load_from_read(reader)
            .expect("Toy ontology should be OK");

        // GENE 1
        let symbol1 = String::from("gene1");

        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0052693".parse().unwrap();

        let mut gene1_hashset = HashSet::new();
        gene1_hashset.insert(t1.clone());
        gene1_hashset.insert(t2.clone());

        let mut term_map = HashMap::new();
        term_map.insert(symbol1.clone(), gene1_hashset);

        // TO DO: add something for tissue expression
        let mut tissue_map = HashMap::new();
        tissue_map.insert(symbol1.clone(), HashSet::new());

        //TO DO: ass something for phenotypes
        let mut phenotypes = HashMap::new();
        phenotypes.insert(symbol1.clone(), HashSet::new());

        let gene_set = Box::leak(Box::new(GeneSetAnnotations::from(term_map, tissue_map, phenotypes)));
        
        NaiveSatisfactionChecker::new(go, gene_set)
    }

    #[test]
    fn test_is_satisfied_esssential() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0052693".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec: Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, false));
        term_vec.push(TermObservation::new(t2, false));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction {
            term_observations: term_vec,
            tissue_expressions: vec![], //to change, in order to test also tissue_expressions
        };

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, true);

        //TO DO: test also tissue expression part
    }

    #[test]
    fn test_is_satisfied_subset() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec: Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, false));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction {
            term_observations: term_vec,
            tissue_expressions: vec![], //to change, in order to test also tissue_expressions
        };

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, true);
    }

    #[test]
    fn test_is_not_satisfied_essential() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0052693".parse().unwrap();
        let t3: TermId = "GO:0005634".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec: Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, false));
        term_vec.push(TermObservation::new(t2, false));
        term_vec.push(TermObservation::new(t3, false));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction {
            term_observations: term_vec,
            tissue_expressions: vec![], //to change, in order to test also tissue_expressions
        };

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, false);
    }

    #[test]
    fn test_is_satisfied_with_exclusion() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0052693".parse().unwrap();
        let t3: TermId = "GO:0005634".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec: Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, false));
        term_vec.push(TermObservation::new(t2, false));
        term_vec.push(TermObservation::new(t3, true));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction {
            term_observations: term_vec,
            tissue_expressions: vec![], //to change, in order to test also tissue_expressions
        };

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, true);
    }

    #[test]
    fn test_is_not_satisfied_with_exclusion() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0052693".parse().unwrap();
        let t3: TermId = "GO:0005634".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec: Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, true));
        term_vec.push(TermObservation::new(t2, false));
        term_vec.push(TermObservation::new(t3, true));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction {
            term_observations: term_vec,
            tissue_expressions: vec![], //to change, in order to test also tissue_expressions
        };

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, false);
    }

    #[test]
    fn test_is_satisfied_with_ontology() {
        // gene1 is annotated with GO_0051146
        // GO:0042692 is a parent of GO_0051146
        // The conjunction containing GO:0042692 should return true because GO_0051146 inherits the parent GO:0042692
        let t1: TermId = "GO:0042692".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec: Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, false));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction {
            term_observations: term_vec,
            tissue_expressions: vec![], //to change, in order to test also tissue_expressions
        };

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, true);
    }



    #[test]
    fn test_is_not_satisfied_with_ontology() {
        // GO:0055007 is a child GO_0051146
        let t1: TermId = "GO:0055007".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec: Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, false));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction {
            term_observations: term_vec,
            tissue_expressions: vec![], //to change, in order to test also tissue_expressions
        };

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, false);
    }

    #[test]
    fn test_is_satisfied_with_ontology_and_exclusion() {
        let t1: TermId = "GO:0055007".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec: Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, true));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction {
            term_observations: term_vec,
            tissue_expressions: vec![], //to change, in order to test also tissue_expressions
        };

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, true);
    }


    #[test]
    fn test_is_not_satisfied_with_ontology_and_exclusion() {
        let t1: TermId = "GO:0042692".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec: Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, true));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction {
            term_observations: term_vec,
            tissue_expressions: vec![], //to change, in order to test also tissue_expressions
        };

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, false);
    }

}
