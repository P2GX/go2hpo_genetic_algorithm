use std::{fs::File, io::BufReader};

use flate2::bufread::GzDecoder;
use ontolius::{
    io::OntologyLoaderBuilder,
    ontology::csr::MinimalCsrOntology, prelude::TermAware,
};

#[test]
fn test_loading_go() {
    let go_path = "data/go/go.toy.json.gz";
    let reader = GzDecoder::new(BufReader::new(
        File::open(go_path).expect("The file should be in the repo"),
    ));

    let parser = OntologyLoaderBuilder::new().obographs_parser().build();
    let go: MinimalCsrOntology = parser
        .load_from_read(reader)
        .expect("The ontology file should be OK");

    assert_eq!(go.len(), 33); // this includes the terms that are not GO terms, e.g. `OIO#hasRelatedSynonym`
}
