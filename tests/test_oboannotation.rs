use clap::Parser;
use std::path::PathBuf;

use ontolius::ontology::csr::MinimalCsrOntology;
use go2hpo_genetic_algorithm::logical_formula::NaiveSatisfactionChecker;

use oboannotation::goannotation::GoAnnotations;


#[test]
fn test_oboannotation () {
    let goa_path_str: &str = "data/gaf/goa_human.gaf.gz";
    println!("processing {}", goa_path_str);
    let goannots = GoAnnotations::new(goa_path_str).expect("Could not create GO Annotations");
    let result = goannots.get_annotation_statistics_json();
    match result {
        Ok(json_string) => println!("{}", json_string),
        Err(e) => eprint!("Could not extract GOA stats: {}", e.to_string())
    }
    // now count the annotations for the first ten genes
    println!("##############\nShowing the first ten GO annotation profiles");
    let annot_map = goannots.get_annotation_map();
    let mut c = 0;
    for (symbol, term_ids) in &annot_map {
        let count = term_ids.len();
        println!("{} has {} unique GO annotations", symbol, count);
        if c > 10 {
            break;
        } else {
            c += 1;
        }
    }
}


#[test]
fn test_gene_annotations () {
    let goa_path_str: &str = "data/gaf/goa_human.gaf.gz";
    println!("processing {}", goa_path_str);
    let goannots = GoAnnotations::new(goa_path_str).expect("Could not create GO Annotations");

    // now count the annotations for the first ten genes
    let annot_map = goannots.get_annotation_map();
    let mut c = 0;
    for (symbol, term_ids) in &annot_map {
        let count = term_ids.len();
        println!("{} has {} unique GO annotations", symbol, count);
        if c > 10 {
            break;
        } else {
            c += 1;
        }
    }

    let checker : NaiveSatisfactionChecker<MinimalCsrOntology>;


}