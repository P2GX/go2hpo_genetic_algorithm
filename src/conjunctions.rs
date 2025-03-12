use std::{cell::OnceCell, collections::{HashMap, HashSet}};
use lazy_static::lazy_static;

use ontolius::{
    base::TermId,
    ontology::csr::MinimalCsrOntology,
    prelude::{
        AncestorNodes, DescendantNodes, HierarchyAware,TermAware,
    },
};

pub struct TermObservation {
    term_id: TermId,
    is_excluded: bool,
}

impl TermObservation{
    pub fn new(term_id: TermId, is_excluded: bool) -> Self {
        Self { term_id, is_excluded }
    }
}

pub struct Conjunction {
    term_observations: Vec<TermObservation>,
}

pub trait SatisfactionChecker {
    fn is_satisfied(&self, symbol: &str, conjunction: &Conjunction) -> bool;
}

pub trait HierarchyTraversal {
    fn get_ancestors(&self, query: &TermId) -> HashSet<TermId>;
    fn get_descendants(&self, query: &TermId) -> HashSet<TermId>;
}

impl HierarchyTraversal for MinimalCsrOntology {
    fn get_ancestors(&self, query: &TermId) -> HashSet<TermId> {
        if let Some(idx) = self.id_to_idx(query) {
            self.hierarchy()
                .iter_ancestors_of(idx)
                .map(|idx| {
                    self.idx_to_term_id(idx).expect(
                        "Ontology should contain a term ID for an index that it just gave up",
                    )
                })
                .cloned()
                .collect()
        } else {
            panic!("query {query} is not in ontology")
        }
    }

    fn get_descendants(&self, query: &TermId) -> HashSet<TermId> {
        if let Some(idx) = self.id_to_idx(query) {
            self.hierarchy()
                .iter_descendants_of(idx)
                .map(|idx| {
                    self.idx_to_term_id(idx).expect(
                        "Ontology should contain a term ID for an index that it just gave up",
                    )
                })
                .cloned()
                .collect()
        } else {
            panic!("query {query} is not in ontology")
        }
    }
}

pub struct NaiveSatisfactionChecker<O> {
    go: O,
    symbol_to_direct_annotations: HashMap<String, HashSet<TermId>>,
}
impl<O> SatisfactionChecker for NaiveSatisfactionChecker<O>
where
    O: HierarchyTraversal,
{
    fn is_satisfied(&self, symbol: &str, conjunction: &Conjunction) -> bool {
        let result = self.symbol_to_direct_annotations.get(symbol);
        match result {
            Some(go_annots_set) => {
                for term_ob in &conjunction.term_observations {

                    // let term_desc: HashSet<TermId> = self.go.get_descendants(&term_ob.term_id);
                    let term_desc: OnceCell<HashSet<TermId>> = OnceCell::new();
                    //LONGER VERSION:
                    match term_ob.is_excluded {
                        true => {
                            if go_annots_set.contains(&term_ob.term_id) {
                                return false;
                            }
                            // else if !go_annots_set.is_disjoint(&term_desc){
                            //     return false;
                            // }

                            // // ONCE CELL ALTERNATIVE
                            let descendants = term_desc.get_or_init(|| self.go.get_descendants(&term_ob.term_id));
                            if !go_annots_set.is_disjoint(descendants){
                                return false;
                            }
                        } 
                        false => {
                            if !go_annots_set.contains(&term_ob.term_id) {
                                let descendants = term_desc.get_or_init(|| self.go.get_descendants(&term_ob.term_id));
                                if go_annots_set.is_disjoint(descendants){
                                    return false;
                                }
                            }
                        }
                    }

                    // // SHORTER VERSION: (PREFERRED)
                    // if term_ob.is_excluded == go_annots_set.contains(&term_ob.term_id){
                    //     return false;
                    // }
                    // else if term_ob.is_excluded == !go_annots_set.is_disjoint(&term_desc){
                    //     return false;
                    // }

                    // // EVEN SHORTER VERSION:
                    // if term_ob.is_excluded == (go_annots_set.contains(&term_ob.term_id) || !go_annots_set.is_disjoint(&term_desc)){
                    //     return false;
                    // }
                }
                true
            }
            None => panic!("We could not find gene symbol: {}.", symbol),
        }
    }
}

impl<O> NaiveSatisfactionChecker<O> {
    pub fn new(go: O, map: HashMap<String, HashSet<TermId>>) -> Self {
        Self {
            go,
            symbol_to_direct_annotations: map,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::{fs::File, hash::Hash, io::BufReader};

    use flate2::bufread::GzDecoder;
    use ontolius::io::OntologyLoaderBuilder;

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
    
        
        let mut map = HashMap::new();
        map.insert(symbol1.clone(), gene1_hashset);
    
        NaiveSatisfactionChecker::new(go, map)
    }


    #[test]
    fn test_is_satisfied_esssential() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0052693".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec:Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, false));
        term_vec.push(TermObservation::new(t2, false));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction{term_observations: term_vec};

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, true);
    }

    #[test]
    fn test_is_satisfied_subset() {
        let t1: TermId = "GO:0051146".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec:Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, false));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction{term_observations: term_vec};

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, true);
    }

    #[test]
    fn test_is_not_satisfied_essential(){
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0052693".parse().unwrap();
        let t3: TermId = "GO:0005634".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec:Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, false));
        term_vec.push(TermObservation::new(t2, false));
        term_vec.push(TermObservation::new(t3, false));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction{term_observations: term_vec};

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, false);
    }

    #[test]
    fn test_is_satisfied_with_exclusion(){
        let t1: TermId = "GO:0051146".parse().unwrap();
        let t2: TermId = "GO:0052693".parse().unwrap();
        let t3: TermId = "GO:0005634".parse().unwrap();
        let symbol = String::from("gene1");

        let mut term_vec:Vec<TermObservation> = Vec::new();
        term_vec.push(TermObservation::new(t1, false));
        term_vec.push(TermObservation::new(t2, false));
        term_vec.push(TermObservation::new(t3, true));
        // term_vec.push(TermObservation::new(t1, false));

        let conjunction = Conjunction{term_observations: term_vec};

        let actual = checker.is_satisfied(&symbol, &conjunction);

        assert_eq!(actual, true);
    }

    // #[test]
    fn test_is_satisfied_with_ontology(){
        
        
    }


    // #[test]
    fn test_is_not_satisfied_with_ontology(){
        todo!()
    }
}
