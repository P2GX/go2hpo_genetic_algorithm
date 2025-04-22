use std::collections::{HashMap, HashSet};
use crate::logical_formula::{DgeState, TissueExpression};
use ontolius::TermId;

use crate::logical_formula::Conjunction;


pub struct GeneAnnotations{
    id: String,
    symbol: String,
    term_annotations: HashSet<TermId>,
    tissue_expressions: HashSet<TissueExpression>,
}

impl GeneAnnotations{
    pub fn contains_term_annotation(&self, term_id: TermId) -> bool{
        todo!()
    }

    pub fn contains_tissue_expressions(&self, tissue_expr: &TissueExpression) -> bool{
        match tissue_expr.state{
            DgeState::Up => self.tissue_expressions.contains(tissue_expr), 
            DgeState::Down => self.tissue_expressions.contains(tissue_expr),
            DgeState::Normal => !self.tissue_expressions.contains(&tissue_expr.into_down()) && !self.tissue_expressions.contains(&tissue_expr.into_up()),
        }
    }
}


pub struct GeneSetAnnotations{
    gene_annotations: HashMap<String, GeneAnnotations>
}

impl GeneSetAnnotations{
    pub fn new(symbol_to_direct_annotations: HashMap<String, HashSet<TermId>>, gene_tissue_expressions: HashMap<String, HashSet<TissueExpression>>) -> Self{
         todo!()
    }

    pub fn get(&self, gene: &str) -> Option<&GeneAnnotations>{
        self.gene_annotations.get(gene)
    }


    /// to check how many keys between the two hashmaps are shared. Useful to detect diffent ID types/formats
    pub fn get_keys_intersection<T, U>(
        set1: &HashMap<String, T>,
        set2: &HashMap<String, U>,
    ) -> HashSet<String> {
        let keys1: HashSet<_> = set1.keys().cloned().collect();
        let keys2: HashSet<_> = set2.keys().cloned().collect();
    
        keys1
            .intersection(&keys2)
            .cloned()
            .collect()
    }

    pub fn get_keys_intersection_percentage<T, U>(
        set1: &HashMap<String, T>,
        set2: &HashMap<String, U>,
    )-> f64 {
        let intersect_size = GeneSetAnnotations::get_keys_intersection(set1, set2).len();
        let smaller_set_size = min(set1.len(), set2.len());
        return intersect_size as f64 / smaller_set_size as f64;
    }
}

