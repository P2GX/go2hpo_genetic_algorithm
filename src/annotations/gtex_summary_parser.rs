use std::collections::{HashMap, HashSet};

use gtex_analyzer::expression_analysis::GtexSummary;

use crate::logical_formula::{DgeState, TissueExpression};

pub struct GtexSummaryParser {}

impl GtexSummaryParser {
    pub fn parse(summary: &GtexSummary) -> HashMap<String, HashSet<TissueExpression>> {
        let mut new_hashmap: HashMap<String, HashSet<TissueExpression>> = HashMap::new();
        let results = summary.get_results();
        for (gene, result) in results {
            let mut tissue_expr_set: HashSet<TissueExpression> = result
                .up_regulated
                .iter()
                .map(|tissue_analysis| {
                    TissueExpression::new(tissue_analysis.tissue_name.clone(), DgeState::Up)
                })
                .collect();

            tissue_expr_set.extend(result.down_regulated.iter().map(|tissue_analysis| {
                TissueExpression::new(tissue_analysis.tissue_name.clone(), DgeState::Down)
            }));

            new_hashmap.insert(gene.to_string(), tissue_expr_set);
        }
        new_hashmap
    }
}

//TO DO: a test for it

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::GzDecoder;
    use gtex_analyzer::expression_analysis::GtexSummaryLoader;
    use std::fs::File;
    use std::io::{self, BufRead, BufReader, Cursor, Read};
    use std::path::Path;

    #[test]
    fn test_gtex_parse() -> io::Result<()> {
        let file_path: &str = "data/gtex/GTEx_RNASeq_gene_median_tpm_HEAD.gct";

        let file = File::open(file_path)?;
        let reader = BufReader::new(file);

        let summary_loader = GtexSummaryLoader::new(Some(10), None);
        let summary = summary_loader.load_summary(reader)?;

        let results = summary.get_results();

        let tissue_expressions = GtexSummaryParser::parse(&summary);

        let mut i = 0;
        for (gene, result) in results {
            //check that it exists
            assert!(tissue_expressions.contains_key(gene));

            let tissue_expr_set = tissue_expressions.get(gene).unwrap();
            if result.up_regulated.len() > 0 {
                let tissue_up = &result.up_regulated.get(0).unwrap().tissue_name;
                let tissue_expr_up = TissueExpression::new(tissue_up.to_string(), DgeState::Up);
                assert!(tissue_expr_set.contains(&tissue_expr_up));
            }

            if result.down_regulated.len() > 0 {
                let tissue_down = &result.down_regulated.get(0).unwrap().tissue_name;
                let tissue_expr_down =
                    TissueExpression::new(tissue_down.to_string(), DgeState::Down);
                assert!(tissue_expr_set.contains(&tissue_expr_down));
            }

            i += 1;
            if i >= 10 {
                break;
            }
        }

        Ok(())
    }
}
