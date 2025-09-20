use flate2::read::GzDecoder;
use serde_json::Value;
use std::fs::File;
use std::io::{BufReader, Read};

pub fn preview_go_ids(path: &str, limit: usize) -> std::io::Result<()> {
    let file = File::open(path)?;
    let mut gz = GzDecoder::new(BufReader::new(file));
    let mut contents = String::new();
    gz.read_to_string(&mut contents)?;

    let json: Value = serde_json::from_str(&contents).expect("Invalid JSON");
    if let Some(graphs) = json.get("graphs").and_then(|g| g.as_array()) {
        for graph in graphs {
            if let Some(nodes) = graph.get("nodes").and_then(|n| n.as_array()) {
                println!("Previewing first {} node IDs:", limit);
                for node in nodes.iter().take(limit) {
                    if let Some(id) = node.get("id").and_then(|id| id.as_str()) {
                        println!("{}", id);
                    }
                }
            }
        }
    }

    Ok(())
}

fn main() -> std::io::Result<()> {
    preview_go_ids("data/go/go-basic.json.gz", 20)?;
    Ok(())
}
