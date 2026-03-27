// src/ontology_io.rs

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};

use anyhow::{Error, Result};
use bincode;
use ontolius::term::AltTermIdAware;
use serde::{Deserialize, Serialize};

use ontolius::io::{GraphEdge as OntGraphEdge, OntologyData, Relationship as OntRelationship};
use ontolius::ontology::csr::MinimalCsrOntology;
use ontolius::term::simple::SimpleMinimalTerm;
use ontolius::term::MinimalTerm;
use ontolius::{Identified, TermId};

//
// ---------- Newtype wrapper to bypass orphan rule ----------
//
pub struct RawOntologyData<I, T>(pub OntologyData<I, T>);

impl<I, T> TryFrom<OntologyData<I, T>> for RawOntologyData<I, T> {
    type Error = Error;
    fn try_from(value: OntologyData<I, T>) -> Result<Self, Self::Error> {
        Ok(RawOntologyData(value))
    }
}

//
// ---------- Serializable wrappers ----------
//

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SerRelationship {
    Parent,
    Child,
    PartOf,
}

impl From<OntRelationship> for SerRelationship {
    fn from(r: OntRelationship) -> Self {
        match r {
            OntRelationship::Parent => SerRelationship::Parent,
            OntRelationship::Child => SerRelationship::Child,
            OntRelationship::PartOf => SerRelationship::PartOf,
        }
    }
}
impl From<SerRelationship> for OntRelationship {
    fn from(r: SerRelationship) -> Self {
        match r {
            SerRelationship::Parent => OntRelationship::Parent,
            SerRelationship::Child => OntRelationship::Child,
            SerRelationship::PartOf => OntRelationship::PartOf,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerGraphEdge<I> {
    pub sub: I,
    pub pred: SerRelationship,
    pub obj: I,
}

impl<I: Clone> From<OntGraphEdge<I>> for SerGraphEdge<I> {
    fn from(e: OntGraphEdge<I>) -> Self {
        Self {
            sub: e.sub,
            pred: e.pred.into(),
            obj: e.obj,
        }
    }
}
impl<I: Clone> From<SerGraphEdge<I>> for OntGraphEdge<I> {
    fn from(e: SerGraphEdge<I>) -> Self {
        Self {
            sub: e.sub,
            pred: e.pred.into(),
            obj: e.obj,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerSimpleMinimalTerm {
    pub term_id: String,
    pub alt_term_ids: Vec<String>,
    pub name: String,
    pub is_obsolete: bool,
}

impl From<SimpleMinimalTerm> for SerSimpleMinimalTerm {
    fn from(term: SimpleMinimalTerm) -> Self {
        SerSimpleMinimalTerm {
            term_id: term.identifier().to_string(), // via Identified
            alt_term_ids: term.iter_alt_term_ids().map(|id| id.to_string()).collect(),
            name: term.name().to_string(),   // via MinimalTerm
            is_obsolete: term.is_obsolete(), // via MinimalTerm
        }
    }
}

impl From<SerSimpleMinimalTerm> for SimpleMinimalTerm {
    fn from(term: SerSimpleMinimalTerm) -> Self {
        let term_id: TermId = term.term_id.parse().expect("valid TermId string");
        let alt_ids: Vec<TermId> = term
            .alt_term_ids
            .into_iter()
            .map(|s| s.parse().expect("valid TermId string"))
            .collect();

        SimpleMinimalTerm::new(term_id, term.name, alt_ids, term.is_obsolete)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerOntologyData<I> {
    pub terms: Vec<SerSimpleMinimalTerm>,
    pub edges: Vec<SerGraphEdge<I>>,
    pub metadata: HashMap<String, String>,
}

impl<'a, I: Clone> From<&'a OntologyData<I, SimpleMinimalTerm>> for SerOntologyData<I> {
    fn from(data: &'a OntologyData<I, SimpleMinimalTerm>) -> Self {
        SerOntologyData {
            // SimpleMinimalTerm is Clone; convert each cloned term to the serializable wrapper
            terms: data.terms.iter().cloned().map(Into::into).collect(),
            // GraphEdge<I> is Clone; convert each to SerGraphEdge<I>
            edges: data.edges.iter().cloned().map(Into::into).collect(),
            metadata: data.metadata.clone(),
        }
    }
}

impl<I: Clone> From<SerOntologyData<I>> for OntologyData<I, SimpleMinimalTerm> {
    fn from(data: SerOntologyData<I>) -> Self {
        Self {
            terms: data.terms.into_iter().map(Into::into).collect(),
            edges: data.edges.into_iter().map(Into::into).collect(),
            metadata: data.metadata,
        }
    }
}

//
// ---------- Save / Load helpers ----------
//

/// Save OntologyData (lossless) via SerOntologyData.
pub fn save_ontology_data<I>(data: &OntologyData<I, SimpleMinimalTerm>, path: &str) -> Result<()>
where
    I: Clone + Serialize + for<'de> Deserialize<'de>,
{
    let ser: SerOntologyData<I> = data.into();
    let file = File::create(path)?;
    let writer = BufWriter::new(file);
    bincode::serialize_into(writer, &ser)?;
    Ok(())
}

/// Load OntologyData back from cache.
pub fn load_ontology_data<I>(path: &str) -> Result<OntologyData<I, SimpleMinimalTerm>>
where
    I: Clone + Serialize + for<'de> Deserialize<'de>,
{
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let ser: SerOntologyData<I> = bincode::deserialize_from(reader)?;
    Ok(ser.into())
}

/// Load MinimalCsrOntology from cached OntologyData.
pub fn load_minimal_csr(path: &str) -> Result<MinimalCsrOntology> {
    let data: OntologyData<u32, SimpleMinimalTerm> = load_ontology_data(path)?;
    let csr: MinimalCsrOntology = data.try_into()?; // works because TryFrom is defined
    Ok(csr)
}
