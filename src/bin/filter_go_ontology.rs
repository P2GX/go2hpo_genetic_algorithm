use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use serde_json::Value;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read};

/// Return true if the ID is a GO term (either starts with "GO:" or contains "GO_").
fn is_go_id(raw: &str) -> bool {
    raw.starts_with("GO:") || raw.contains("GO_")
}

/// Filter ontology JSON to keep only GO terms.
/// Keeps only nodes/edges whose IDs look like GO terms.
pub fn filter_go_json(input_path: &str, output_path: &str) -> std::io::Result<()> {
    // Open gzipped input
    let file = File::open(input_path)?;
    let mut gz = GzDecoder::new(BufReader::new(file));
    let mut contents = String::new();
    gz.read_to_string(&mut contents)?;

    // Parse JSON
    let mut json: Value = serde_json::from_str(&contents)?;

    if let Some(graphs) = json.get_mut("graphs").and_then(|g| g.as_array_mut()) {
        for graph in graphs {
            // --- Filter nodes ---
            if let Some(nodes) = graph.get_mut("nodes").and_then(|n| n.as_array_mut()) {
                nodes.retain(|node| {
                    node.get("id")
                        .and_then(|id| id.as_str())
                        .map(is_go_id)
                        .unwrap_or(false)
                });
            }

            // --- Filter edges ---
            if let Some(edges) = graph.get_mut("edges").and_then(|e| e.as_array_mut()) {
                edges.retain(|edge| {
                    let subj_ok = edge
                        .get("sub")
                        .and_then(|id| id.as_str())
                        .map(is_go_id)
                        .unwrap_or(false);
                    let obj_ok = edge
                        .get("obj")
                        .and_then(|id| id.as_str())
                        .map(is_go_id)
                        .unwrap_or(false);
                    subj_ok && obj_ok
                });
            }
        }
    }

    // Save gzipped output
    let out_file = File::create(output_path)?;
    let gz_out = GzEncoder::new(BufWriter::new(out_file), Compression::default());
    serde_json::to_writer(gz_out, &json)?;

    Ok(())
}

fn main() -> std::io::Result<()> {
    filter_go_json(
        "data/go/go-basic.json.gz",
        "data/go/go-basic.filtered.json.gz",
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::bufread::GzDecoder;
    use ontolius::io::OntologyLoaderBuilder;
    use ontolius::ontology::csr::MinimalCsrOntology;
    use ontolius::ontology::OntologyTerms;
    use std::fs::File;
    use std::io::BufReader;

    fn load_ontology(path: &str) -> MinimalCsrOntology {
        let reader = GzDecoder::new(BufReader::new(
            File::open(path).expect("The file should be available. File not found."),
        ));

        let parser = OntologyLoaderBuilder::new().obographs_parser().build();
        parser
            .load_from_read(reader)
            .expect("The ontology file should be OK")
    }

    #[test]
    fn go_filtered_vs_unfiltered_counts() {
        let unfiltered_path = "data/go/go-basic.json.gz";
        let filtered_path = "data/go/go-basic.filtered.json.gz";

        let go_unfiltered = load_ontology(unfiltered_path);
        let go_filtered = load_ontology(filtered_path);

        let terms_unfiltered = go_unfiltered.iter_term_ids().count();
        let terms_filtered = go_filtered.iter_term_ids().count();
        let all_terms_unfiltered = go_unfiltered.iter_all_term_ids().count();
        let all_terms_filtered = go_filtered.iter_all_term_ids().count();

        println!("Unfiltered ontology:");
        println!("  iter_term_ids():     {}", terms_unfiltered);
        println!("  iter_all_term_ids(): {}", all_terms_unfiltered);

        println!("Filtered ontology:");
        println!("  iter_term_ids():     {}", terms_filtered);
        println!("  iter_all_term_ids(): {}", all_terms_filtered);

        // Just a sanity check: filtered ontology must have less or equal terms
        assert!(terms_filtered <= terms_unfiltered);
        assert!(all_terms_filtered <= all_terms_unfiltered);
    }
}
