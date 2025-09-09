use anyhow::bail;
use flate2::read::GzDecoder;
use ontolius::io::OntologyLoaderBuilder;
use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;
use std::{fs::File, io::BufReader};

use go2hpo_genetic_algorithm::logical_formula::NaiveSatisfactionChecker;
use go2hpo_genetic_algorithm::logical_formula::SatisfactionChecker;
use go2hpo_genetic_algorithm::annotations::GeneSetAnnotations;
use ontolius::ontology::csr::MinimalCsrOntology;

use oboannotation::go::GoAnnotations;

use oboannotation::go::GoGafAnnotationLoader;
use oboannotation::io::AnnotationLoader;

use oboannotation::go::stats::get_annotation_map;

use go2hpo_genetic_algorithm::logical_formula::Conjunction;

use go2hpo_genetic_algorithm::logical_formula::TermObservation;

use ontolius::TermId;
use crate::fixtures::gene_set_annotations::gene_set_annotations;
use lazy_static::lazy_static;
use rstest::{fixture, rstest};
mod fixtures;

lazy_static! {
    static ref ANNOTATION_MAP: HashMap<String, HashSet<TermId>>  = load_and_get_go_annotation_map();
}



#[test]
fn test_oboannotation() -> anyhow::Result<()> {
    // now count the annotations for the first ten genes
    println!("##############\nShowing the first ten GO annotation profiles");
    let annotations = load_and_get_go_annotations().expect("Load GoAnnotations");
    let annot_map = get_annotation_map(&annotations);
    let mut c = 0;
    for (symbol, term_ids) in &annot_map {
        println!("{} has {} unique GO annotations", symbol, term_ids.len());
        if c > 10 {
            break;
        } else {
            c += 1;
        }
    }
    Ok(())
}

fn load_and_get_go_annotation_map() -> HashMap<String, HashSet<TermId>>{
    let annotations = load_and_get_go_annotations().expect("Load GoAnnotations");
    get_annotation_map(&annotations)
}
fn load_and_get_go_annotations() -> anyhow::Result<GoAnnotations>{
    let goa_path_str: &str = "data/gaf/goa_human.gaf.gz";
    println!("processing {}", goa_path_str);
    let reader: Box<dyn BufRead> = open_for_reading(goa_path_str)?;

    let loader = GoGafAnnotationLoader;
    let annotations = match loader.load_from_buf_read(reader) {
        Ok(annotations) => {
            println!(
                "Loaded {:?} annotations. Skipped {:?} negated annotations",
                annotations.annotations.len(),
                annotations.negated_annotation_count
            );
            annotations
        }
        Err(e) => {
            bail!("Could not load GOA: {e}")
        }
    };
    Ok(annotations)
}


#[test]
fn show_some(){
    let annot_map = ANNOTATION_MAP.clone();
    let mut c = 0;
    for (symbol, term_ids) in &annot_map {
        if term_ids.len() > 2 && term_ids.len() < 10{
            println!("{} has {} unique GO annotations", symbol, term_ids.len());
            println!("\t {}: {:?}", symbol, term_ids);
            if c > 10 {
                break;
            } else {
                c += 1;
            }
        }
    }
}

#[test]
fn show_pmpca_gene(){
    let annot_map = ANNOTATION_MAP.clone();
    let symbol = String::from("PMPCA");
    let term_ids = annot_map.get(&symbol).unwrap();
    println!("{}: {:?}", symbol, term_ids);
}


fn load_and_get_go() -> MinimalCsrOntology{
    let go_path = "data/go/go.toy.json.gz";
    let reader = GzDecoder::new(BufReader::new(
        File::open(go_path).expect("The file should be in the repo"),
    ));
    let loader = OntologyLoaderBuilder::new().obographs_parser().build();
    let go: MinimalCsrOntology = loader
        .load_from_read(reader)
        .expect("Toy ontology should be OK");
    go
}

macro_rules! test_conjunction {
    ($func:ident, $symbol:expr, [$(($go:expr, $val:expr)),+], $expected:expr) => {
        
        #[rstest]
        fn $func(gene_set_annotations: GeneSetAnnotations){
            let go: MinimalCsrOntology = load_and_get_go();
            let annot_map = ANNOTATION_MAP.clone();
            let gene_set_annotations: &'static GeneSetAnnotations = Box::leak(Box::new(gene_set_annotations));
            let checker: NaiveSatisfactionChecker<MinimalCsrOntology> = NaiveSatisfactionChecker::new(&go, &gene_set_annotations);

            let mut term_vec: Vec<TermObservation> = Vec::new();
            $(
                let term: TermId = $go.parse().unwrap();
                term_vec.push(TermObservation::new(term, $val));
            )+

            let conjunction = Conjunction {
                term_observations: term_vec,
                tissue_expressions: vec![],
            };

            let result = checker.is_satisfied($symbol, &conjunction);
            assert_eq!(result, $expected);
            
            
    }};
}

test_conjunction!(test_conjunction_not_satisfied , &"PMPCA".to_string(), [("GO:0051146", false), ("GO:0052693", false)], false);
test_conjunction!(test_conjunction_satisfied,  &"PMPCA".to_string(), [("GO:0051146", true), ("GO:0052693", true)], true);
test_conjunction!(test_conjunction_satisfied_with_included_terms,  &"PMPCA".to_string(), [("GO:0004222", false), ("GO:0017087", false)], true);


fn open_for_reading<P: AsRef<Path>>(goa_path: P) -> anyhow::Result<Box<dyn BufRead>> {
    Ok(if let Some(extension) = goa_path.as_ref().extension() {
        if extension == "gz" {   
            Box::new(BufReader::new(GzDecoder::new(File::open(goa_path)?)))
        } else {
            Box::new(BufReader::new(File::open(goa_path)?))
        }
    } else {
        Box::new(BufReader::new(File::open(goa_path)?))
    })
}
