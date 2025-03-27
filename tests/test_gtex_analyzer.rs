use flate2::read::GzDecoder;
use gtex_analyzer::expression_analysis::GtexSummaryLoader;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::path::Path;



#[test]
fn test_load_gtex_analyzer() -> io::Result<()> {
    let file_path: &str = "data/gtex/GTEx_RNASeq_gene_median_tpm_HEAD.gct";

    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    let summary_loader = GtexSummaryLoader::new(Some(10), None);
    let summary = summary_loader.load_summary(reader)?;

    // println!("{:#?}", summary.get_results());

    // assert!(summary.metadata.is_some(), "Metadata should be present");
    assert!(
        !summary.get_results().is_empty(),
        "Results should not be empty"
    );

    assert!(summary.metadata.num_tissues > 0);

    Ok(())
}