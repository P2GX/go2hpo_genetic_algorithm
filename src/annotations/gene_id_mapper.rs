use std::collections::HashMap;

use super::GeneId;


pub struct GeneIdMapper{
    map: HashMap<String, GeneId>,
}

impl GeneIdMapper{
    pub fn new(map: HashMap<String, GeneId>) -> Self{
        Self { map }
    }
    pub fn map(&self, term: &str) -> Option<GeneId>{
        self.map.get(term).cloned()
    }

    pub fn map_keys<T>(&self, hashmap: HashMap<String, T>) -> HashMap<GeneId, T>{
        let mut new_hashmap: HashMap<GeneId, T> = HashMap::new();
        let mut n_unmapped = 0;
        for (old_gene_id, content) in hashmap{
            match self.map(&old_gene_id){
                Some(new_gene_id) => {new_hashmap.insert(new_gene_id, content);},
                None => {n_unmapped += 1;},
            }
        }

        if n_unmapped > 0 {println!("A total of {} terms haven't been mapped successfully", n_unmapped)}

        new_hashmap
    }
}

#[cfg(test)]
mod tests {
    use gtex_analyzer::expression_analysis::GtexSummaryLoader;

    use crate::annotations::GtexSummaryParser;

    use super::*;
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{self, BufRead, BufReader};

    fn get_map() -> io::Result<HashMap<String, GeneId>>{
        let file = File::open("data/map/hgnc_complete_set.txt")?;
        let reader = BufReader::new(file);

        let mut map: HashMap<String, GeneId> = HashMap::new();

        let mut lines = reader.lines();

        let header = lines.next().unwrap()?; 
        let columns: Vec<&str> = header.split('\t').collect();

        let ensembl_index = columns.iter().position(|&col| col == "ensembl_gene_id").unwrap();
        let symbol_index = columns.iter().position(|&col| col == "symbol").unwrap();

        for line in lines {
            if let Ok(line) = line {
                let fields: Vec<&str> = line.split('\t').collect();
                if let (Some(ensembl), Some(symbol)) = (fields.get(ensembl_index), fields.get(symbol_index)) {
                    if !ensembl.is_empty() && !symbol.is_empty() {
                        map.insert(ensembl.to_string(), symbol.to_string());
                    }
                }
            }
        }
        Ok(map)
    }
    
    #[test]
    fn test_map(){
        let map = get_map().unwrap();

        println!("Loaded {} mappings", map.len());

        assert!(!map.is_empty(), "The map should not be empty");
    }

    fn remove_ensembl_id_version_number<T>(hashmap: HashMap<String, T>) -> HashMap<String, T>{
        let mut new_hashmap: HashMap<String, T> = HashMap::new();
        for (ensembl_id, content) in hashmap{
            let base_id = ensembl_id.split('.').next().unwrap();
            new_hashmap.insert(base_id.to_string(), content);
        }   
        new_hashmap
    }
    
    #[test]
    fn test_gene_id_mapper() -> io::Result<()>{
        let file_path: &str = "data/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz";
        let file = File::open(file_path)?;
        let gz = flate2::read::GzDecoder::new(file);
        let reader = BufReader::new(gz);
    
        let summary_loader = GtexSummaryLoader::new(None, None);
        let summary = summary_loader.load_summary(reader)?;
        
        let tissue_expressions = remove_ensembl_id_version_number(GtexSummaryParser::parse(&summary));

        let gene_id_mapper = GeneIdMapper::new(get_map().unwrap());
        
        let mapped_tissue_expressions = gene_id_mapper.map_keys(tissue_expressions);

        println!("Loaded {} mappings", mapped_tissue_expressions.len());

        assert!(!mapped_tissue_expressions.is_empty(), "The map should not be empty");

        Ok(())
    }
}