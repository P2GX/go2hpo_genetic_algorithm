//! Satisfaction checker for pre-expanded GO annotations (direct + ancestors).
//! Assumes that each gene's `term_annotations` already contains all inherited terms,
//! so no ontology traversal or caching is needed at evaluation time.

use std::collections::HashMap;

use crate::annotations::{GeneAnnotations, GeneId, GeneSetAnnotations};

use super::{Conjunction, SatisfactionChecker};

pub struct PreexpandedSatisfactionChecker<'a> {
    gene_set_annotations: &'a GeneSetAnnotations,
}

impl<'a> PreexpandedSatisfactionChecker<'a> {
    pub fn new(gene_set_annotations: &'a GeneSetAnnotations) -> Self {
        Self {
            gene_set_annotations,
        }
    }

    #[inline]
    fn are_term_annotations_satisfied(
        &self,
        gene_annotations: &GeneAnnotations,
        conjunction: &Conjunction,
    ) -> bool {
        let annotations = gene_annotations.get_term_annotations();
        for term_ob in &conjunction.term_observations {
            if term_ob.is_excluded {
                if annotations.contains(&term_ob.term_id) {
                    return false;
                }
            } else if !annotations.contains(&term_ob.term_id) {
                return false;
            }
        }
        true
    }

    #[inline]
    fn are_tissue_expressions_satisfied(
        &self,
        gene_annotations: &GeneAnnotations,
        conjunction: &Conjunction,
    ) -> bool {
        conjunction
            .tissue_expressions
            .iter()
            .all(|tissue_expr| gene_annotations.contains_tissue_expressions(tissue_expr))
    }
}

impl<'a> SatisfactionChecker for PreexpandedSatisfactionChecker<'a> {
    fn is_satisfied(&self, symbol: &GeneId, conjunction: &Conjunction) -> bool {
        match self.gene_set_annotations.get_gene_annotations(symbol) {
            Some(gene_annotations) => {
                self.are_term_annotations_satisfied(gene_annotations, conjunction)
                    && self.are_tissue_expressions_satisfied(gene_annotations, conjunction)
            }
            None => panic!("Couldn't find the gene ID in the Gene Annotation Set"),
        }
    }

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

    fn get_gene_set(&self) -> &GeneSetAnnotations {
        self.gene_set_annotations
    }
}
