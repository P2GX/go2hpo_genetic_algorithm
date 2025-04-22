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

    pub fn contains(gene: String) -> bool{
        todo!()
    }

    pub fn get(gene: String) -> GeneAnnotations{
        todo!()
    }
}

