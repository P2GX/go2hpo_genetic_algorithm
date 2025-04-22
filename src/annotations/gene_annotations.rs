use std::{cmp::min, collections::{HashMap, HashSet}};
use crate::logical_formula::{DgeState, TissueExpression};
use ontolius::TermId;

use crate::logical_formula::Conjunction;


pub struct GeneAnnotations{
    id: String,
    // symbol: String,
    term_annotations: HashSet<TermId>,
    tissue_expressions: HashSet<TissueExpression>,
}

impl GeneAnnotations{
    pub fn new(id: String, term_annotations: HashSet<TermId>, tissue_expressions: HashSet<TissueExpression>) -> Self{
        Self {id, term_annotations, tissue_expressions}
    }
    pub fn contains_term_annotation(&self, term_id: &TermId) -> bool{
        self.term_annotations.contains(term_id)
    }

    pub fn contains_tissue_expressions(&self, tissue_expr: &TissueExpression) -> bool{
        match tissue_expr.state{
            DgeState::Up => self.tissue_expressions.contains(tissue_expr), 
            DgeState::Down => self.tissue_expressions.contains(tissue_expr),
            DgeState::Normal => !self.tissue_expressions.contains(&tissue_expr.into_down()) && !self.tissue_expressions.contains(&tissue_expr.into_up()),
        }
    }

    pub fn get_term_annotations(&self) -> &HashSet<TermId>{
        return &self.term_annotations;
    }

    pub fn get_tissue_expressions(&self) -> &HashSet<TissueExpression>{
        return &self.tissue_expressions;
    }
}


pub struct GeneSetAnnotations{
    gene_annotations: HashMap<String, GeneAnnotations>
}

impl GeneSetAnnotations{
    pub fn from(symbol_to_direct_annotations: HashMap<String, HashSet<TermId>>, gene_tissue_expressions: HashMap<String, HashSet<TissueExpression>>) -> Self{
        // let keys = GeneSetAnnotations::get_keys_intersection(&symbol_to_direct_annotations, &gene_tissue_expressions);
        let mut gene_annotations: HashMap<String, GeneAnnotations> = HashMap::new();
        let intersect_annotations = symbol_to_direct_annotations.iter()
                                .map(|(gene, terms)| (gene, terms, gene_tissue_expressions.get(gene)))
                                .filter(|(gene, terms, tissues)| tissues.is_some())
                                .map(|(gene, terms, tissues)| (gene, terms, tissues.unwrap()));

        for (gene, terms, tissues) in intersect_annotations{
            gene_annotations.insert(gene.to_string(), GeneAnnotations::new(gene.to_string(), terms.clone(), tissues.clone()));
        }

        Self { gene_annotations }
    }

    pub fn contains_key(&self, gene: &str) -> bool{
        self.gene_annotations.contains_key(gene)
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

