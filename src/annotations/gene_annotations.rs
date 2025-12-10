use crate::logical_formula::{DgeState, TermObservation, TissueExpression};
use ontolius::TermId;
use rand::seq::IteratorRandom;
use rand::Rng;
use serde::{Deserialize, Serialize};
use serde_with::{serde_as, DisplayFromStr};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::{
    cmp::min,
    collections::{HashMap, HashSet},
    path::Path,
};

use crate::logical_formula::Conjunction;

pub type GeneId = String;

#[serde_as]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneAnnotations {
    id: GeneId,
    // symbol: String,

    // Annotations
    #[serde_as(as = "HashSet<DisplayFromStr>")]
    term_annotations: HashSet<TermId>,

    tissue_expressions: HashSet<TissueExpression>,

    // Phenotype (HPO)
    #[serde_as(as = "HashSet<DisplayFromStr>")]
    phenotypes: HashSet<TermId>,
}

impl GeneAnnotations {
    pub fn new(
        id: GeneId,
        term_annotations: HashSet<TermId>,
        tissue_expressions: HashSet<TissueExpression>,
        phenotypes: HashSet<TermId>,
    ) -> Self {
        Self {
            id,
            term_annotations,
            tissue_expressions,
            phenotypes,
        }
    }
    pub fn contains_term_annotation(&self, term_id: &TermId) -> bool {
        self.term_annotations.contains(term_id)
    }

    //TO DO: to test
    pub fn contains_tissue_expressions(&self, tissue_expr: &TissueExpression) -> bool {
        match tissue_expr.state {
            DgeState::Up => self.tissue_expressions.contains(tissue_expr),
            DgeState::Down => self.tissue_expressions.contains(tissue_expr),
            DgeState::Normal => {
                !self.tissue_expressions.contains(&tissue_expr.into_down())
                    && !self.tissue_expressions.contains(&tissue_expr.into_up())
            }
        }
    }

    pub fn contains_phenotype(&self, phenotype: &TermId) -> bool {
        self.phenotypes.contains(phenotype)
    }

    pub fn get_term_annotations(&self) -> &HashSet<TermId> {
        return &self.term_annotations;
    }

    pub fn get_tissue_expressions(&self) -> &HashSet<TissueExpression> {
        return &self.tissue_expressions;
    }

    pub fn get_phenotypes(&self) -> &HashSet<TermId> {
        return &self.phenotypes;
    }

    pub fn into_conjunction(&self) -> Conjunction {
        let term_obs: Vec<TermObservation> = self
            .term_annotations
            .iter()
            .map(|term_id| TermObservation::new(term_id.clone(), false))
            .collect();
        let tissue_exprs: Vec<TissueExpression> = self.tissue_expressions.iter().cloned().collect();
        Conjunction::from(term_obs, tissue_exprs)
    }

    // pub fn into_randomly_sampled_conjunction<'a, R: Rng>(&self, prob_terms:f64, prob_tissues:f64, rng: &'a mut R) -> Conjunction{
    //     let term_obs: Vec<TermObservation> = self.term_annotations.iter()
    //                                         .filter(|_| prob_terms > rng.random::<f64>())
    //                                         .map(|term_id| TermObservation::new(term_id.clone(), false))
    //                                         .collect();
    //     let tissue_exprs: Vec<TissueExpression> = self.tissue_expressions.iter().filter(|_| prob_tissues > rng.random::<f64>()).cloned().collect();
    //     Conjunction::from(term_obs, tissue_exprs)
    // }

    pub fn into_randomly_sampled_conjunction<'a, R: Rng>(
        &self,
        prob_terms: f64,
        prob_tissues: f64,
        rng: &'a mut R,
        max_terms: Option<usize>,
        max_tissues: Option<usize>,
    ) -> Conjunction {
        // Step 1: Probabilistic filtering (original behavior)
        let mut term_obs: Vec<TermObservation> = self
            .term_annotations
            .iter()
            .filter(|_| prob_terms > rng.random::<f64>())
            .map(|term_id| TermObservation::new(term_id.clone(), false))
            .collect();

        let mut tissue_exprs: Vec<TissueExpression> = self
            .tissue_expressions
            .iter()
            .filter(|_| prob_tissues > rng.random::<f64>())
            .cloned()
            .collect();

        // Step 2: Apply caps (if provided)
        if let Some(max) = max_terms {
            if term_obs.len() > max {
                term_obs = term_obs
                    .iter()
                    .choose_multiple(rng, max)
                    .into_iter()
                    .cloned()
                    .collect();
            }
        }

        if let Some(max) = max_tissues {
            if tissue_exprs.len() > max {
                tissue_exprs = tissue_exprs
                    .iter()
                    .choose_multiple(rng, max)
                    .into_iter()
                    .cloned()
                    .collect();
            }
        }

        Conjunction::from(term_obs, tissue_exprs)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneSetAnnotations {
    gene_annotations: HashMap<GeneId, GeneAnnotations>,
}

impl GeneSetAnnotations {
    pub fn from(
        symbol_to_direct_annotations: HashMap<GeneId, HashSet<TermId>>,
        gene_tissue_expressions: HashMap<GeneId, HashSet<TissueExpression>>,
        genes_to_phenotype: HashMap<GeneId, HashSet<TermId>>,
    ) -> Self {
        // let keys = GeneSetAnnotations::get_keys_intersection(&symbol_to_direct_annotations, &gene_tissue_expressions);
        let mut gene_annotations: HashMap<String, GeneAnnotations> = HashMap::new();
        let intersect_annotations = symbol_to_direct_annotations
            .iter()
            .map(|(gene, terms)| {
                (
                    gene,
                    terms,
                    gene_tissue_expressions.get(gene),
                    genes_to_phenotype.get(gene),
                )
            })
            .filter(|(gene, terms, tissues, phenotypes)| tissues.is_some() && phenotypes.is_some())
            .map(|(gene, terms, tissues, phenotypes)| {
                (gene, terms, tissues.unwrap(), phenotypes.unwrap())
            });

        for (gene, terms, tissues, phenotypes) in intersect_annotations {
            gene_annotations.insert(
                gene.to_string(),
                GeneAnnotations::new(
                    gene.to_string(),
                    terms.clone(),
                    tissues.clone(),
                    phenotypes.clone(),
                ),
            );
        }

        Self { gene_annotations }
    }

    pub fn new(gene_annotations: HashMap<GeneId, GeneAnnotations>) -> Self {
        Self { gene_annotations }
    }

    pub fn get_gene_annotations_map(&self) -> &HashMap<String, GeneAnnotations> {
        &self.gene_annotations
    }

    // pub fn iter(&self) -> hash_map::Iter<'_, String, GeneAnnotations>{
    //     self.gene_annotations.iter()
    // }

    pub fn len(&self) -> usize {
        return self.gene_annotations.len();
    }

    pub fn contains_gene(&self, gene: &GeneId) -> bool {
        self.gene_annotations.contains_key(gene)
    }

    pub fn get_gene_annotations(&self, gene: &GeneId) -> Option<&GeneAnnotations> {
        self.gene_annotations.get(gene)
    }

    /// Get all GO terms annotated to genes that have the specified HPO phenotype
    pub fn get_go_terms_for_hpo_phenotype(&self, hpo_term: &TermId) -> HashSet<TermId> {
        self.gene_annotations
            .values()
            .filter(|ann| ann.contains_phenotype(hpo_term))
            .flat_map(|ann| ann.get_term_annotations())
            .cloned()
            .collect()
    }

    /// to check how many keys between the two hashmaps are shared. Useful to detect diffent ID types/formats
    pub fn get_keys_intersection<T, U>(
        set1: &HashMap<String, T>,
        set2: &HashMap<String, U>,
    ) -> HashSet<String> {
        let keys1: HashSet<_> = set1.keys().cloned().collect();
        let keys2: HashSet<_> = set2.keys().cloned().collect();

        keys1.intersection(&keys2).cloned().collect()
    }

    pub fn get_keys_intersection_percentage<T, U>(
        set1: &HashMap<String, T>,
        set2: &HashMap<String, U>,
    ) -> f64 {
        let intersect_size = GeneSetAnnotations::get_keys_intersection(set1, set2).len();
        let smaller_set_size = min(set1.len(), set2.len());
        return intersect_size as f64 / smaller_set_size as f64;
    }
}

// GeneSetAnnotations IMPORT / EXPORT

impl GeneSetAnnotations {
    // BINCODE for fast and compact import/export

    pub fn save_bincode<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let file = File::create(path)?;
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
    }

    pub fn load_bincode<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        bincode::deserialize_from(reader)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
    }

    // JSON for human-readable import/export
    pub fn save_json<P: AsRef<Path>>(&self, path: P) -> std::io::Result<()> {
        let file = File::create(path)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, self)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
    }

    pub fn load_json<P: AsRef<Path>>(path: P) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        serde_json::from_reader(reader)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))
    }
}

// TO DO TESTS
