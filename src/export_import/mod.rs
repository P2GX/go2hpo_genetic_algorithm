use std::{collections::HashSet, fs, path::Path, str::FromStr};

use anyhow::{bail, Context, Result};
use serde::{Deserialize, Serialize};

use crate::{
    genetic_algorithm::{FormulaEvaluator, Solution},
    logical_formula::{Conjunction, DNF, DNFVec, TermObservation, TissueExpression},
};
use ontolius::{
    ontology::{csr::MinimalCsrOntology, OntologyTerms},
    TermId,
};

pub const SNAPSHOT_SCHEMA_VERSION: u8 = 1;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableTermObservation {
    pub term_id: String,
    pub is_excluded: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableConjunction {
    pub term_observations: Vec<SerializableTermObservation>,
    pub tissue_expressions: Vec<TissueExpression>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableDNFVec {
    pub conjunctions: Vec<SerializableConjunction>,
}

impl From<&Conjunction> for SerializableConjunction {
    fn from(conj: &Conjunction) -> Self {
        Self {
            term_observations: conj
                .term_observations
                .iter()
                .map(|obs| SerializableTermObservation {
                    term_id: obs.term_id.to_string(),
                    is_excluded: obs.is_excluded,
                })
                .collect(),
            tissue_expressions: conj.tissue_expressions.clone(),
        }
    }
}

impl TryFrom<&SerializableConjunction> for Conjunction {
    type Error = anyhow::Error;

    fn try_from(value: &SerializableConjunction) -> Result<Self> {
        let term_observations = value
            .term_observations
            .iter()
            .map(|obs| {
                let term_id = TermId::from_str(&obs.term_id)
                    .with_context(|| format!("invalid TermId {}", obs.term_id))?;
                Ok(TermObservation::new(term_id, obs.is_excluded))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(Conjunction::from(
            term_observations,
            value.tissue_expressions.clone(),
        ))
    }
}

impl From<&DNFVec> for SerializableDNFVec {
    fn from(dnf: &DNFVec) -> Self {
        Self {
            conjunctions: dnf
                .get_active_conjunctions()
                .iter()
                .map(|c| SerializableConjunction::from(*c))
                .collect(),
        }
    }
}

impl TryFrom<&SerializableDNFVec> for DNFVec {
    type Error = anyhow::Error;

    fn try_from(value: &SerializableDNFVec) -> Result<Self> {
        let conjunctions = value
            .conjunctions
            .iter()
            .map(Conjunction::try_from)
            .collect::<Result<Vec<_>>>()?;
        Ok(DNFVec::from_conjunctions(conjunctions))
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerializableSolution {
    pub formula: SerializableDNFVec,
    pub score: f64,
}

impl From<&Solution<DNFVec>> for SerializableSolution {
    fn from(value: &Solution<DNFVec>) -> Self {
        Self {
            formula: SerializableDNFVec::from(value.get_formula()),
            score: value.get_score(),
        }
    }
}

impl SerializableSolution {
    pub fn to_solution<'a>(
        &self,
        evaluator: &FormulaEvaluator<'a, DNFVec, TermId>,
        phenotype: &TermId,
    ) -> Result<Solution<DNFVec>> {
        let dnf: DNFVec = (&self.formula).try_into()?;
        Ok(evaluator.evaluate(&dnf, phenotype))
    }

    pub fn to_dnf(&self) -> Result<DNFVec> {
        let dnf: DNFVec = (&self.formula).try_into()?;
        Ok(dnf)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaRunMetadata {
    pub schema_version: u8,
    pub hpo_term: String,
    pub pop_size: usize,
    pub generations: usize,
    pub mutation_rate: f64,
    pub tournament_size: usize,
    pub max_n_terms: usize,
    pub max_n_conj: usize,
    pub penalty_lambda: f64,
    pub fscore_beta: f64,
    pub rng_seed: u64,
    pub generation_index: usize,
    pub filtered_go_terms: Option<Vec<String>>,
    pub tissue_terms_used: Vec<String>,
    pub best_score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaPopulationSnapshot {
    pub metadata: GaRunMetadata,
    pub population: Vec<SerializableSolution>,
}

impl GaPopulationSnapshot {
    pub fn best(&self) -> Option<&SerializableSolution> {
        self.population.iter().max_by(|a, b| {
            a.score
                .partial_cmp(&b.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }
}

pub fn export_snapshot<P: AsRef<Path>>(path: P, snapshot: &GaPopulationSnapshot) -> Result<()> {
    let path_ref = path.as_ref();
    if let Some(parent) = path_ref.parent() {
        fs::create_dir_all(parent).context("creating snapshot parent directory")?;
    }

    let bytes = bincode::serialize(snapshot).context("serializing snapshot")?;
    fs::write(path_ref, bytes).context("writing snapshot to disk")?;
    Ok(())
}

pub fn import_snapshot<P: AsRef<Path>>(path: P) -> Result<GaPopulationSnapshot> {
    let data = fs::read(path.as_ref()).context("reading snapshot file")?;
    let snapshot: GaPopulationSnapshot =
        bincode::deserialize(&data).context("deserializing snapshot")?;

    if snapshot.metadata.schema_version != SNAPSHOT_SCHEMA_VERSION {
        bail!(
            "snapshot schema {} is incompatible with reader schema {}",
            snapshot.metadata.schema_version,
            SNAPSHOT_SCHEMA_VERSION
        );
    }

    Ok(snapshot)
}

pub fn rebuild_population_from_snapshot<'a>(
    snapshot: &GaPopulationSnapshot,
    evaluator: &FormulaEvaluator<'a, DNFVec, TermId>,
    phenotype: &TermId,
    go_ontology: &MinimalCsrOntology,
    tissues_available: &[String],
) -> Result<Vec<Solution<DNFVec>>> {
    if snapshot.population.is_empty() {
        bail!("snapshot population is empty");
    }

    let tissue_set: HashSet<&str> = tissues_available.iter().map(|s| s.as_str()).collect();

    let mut population: Vec<Solution<DNFVec>> = Vec::with_capacity(snapshot.population.len());

    for ser_sol in &snapshot.population {
        let dnf = ser_sol.to_dnf()?;

        for conj in dnf.get_active_conjunctions() {
            for obs in &conj.term_observations {
                if go_ontology.term_by_id(&obs.term_id).is_none() {
                    bail!("GO term {} not found in current ontology", obs.term_id);
                }
            }

            for texpr in &conj.tissue_expressions {
                if !tissue_set.contains(texpr.term_id.as_str()) {
                    bail!(
                        "Tissue {} not found in current GTEx metadata",
                        texpr.term_id
                    );
                }
            }
        }

        population.push(evaluator.evaluate(&dnf, phenotype));
    }

    Ok(population)
}
