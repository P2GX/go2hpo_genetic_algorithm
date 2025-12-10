import gzip, csv
from collections import defaultdict
from statistics import mean, median
from pathlib import Path

root = Path('.')
hpo_file = root / 'data/hpo2gene/phenotype_to_genes.txt'
hpo_to_genes = defaultdict(set)
hpo_names = {}
with hpo_file.open('r', encoding='utf-8') as f:
    next(f)
    for line in f:
        hpo_id, hpo_name, ncbi, symbol, disease = line.rstrip('\n').split('\t')
        hpo_to_genes[hpo_id].add(symbol)
        hpo_names[hpo_id] = hpo_name

# GO genes
goa_path = root / 'data/gaf/goa_human.gaf.gz'
go_genes = set()
with gzip.open(goa_path, 'rt', encoding='utf-8', errors='ignore') as f:
    for line in f:
        if line.startswith('!'):
            continue
        parts = line.split('\t')
        if len(parts) > 2:
            symbol = parts[2].strip()
            if symbol:
                go_genes.add(symbol)

# HGNC map
map_path = root / 'data/map/hgnc_complete_set.txt'
ensembl_to_symbol = {}
with map_path.open('r', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter='\t')
    header = next(reader)
    col_idx = {name: i for i, name in enumerate(header)}
    ens_idx = col_idx.get('ensembl_gene_id')
    sym_idx = col_idx.get('symbol')
    for row in reader:
        if len(row) <= max(ens_idx, sym_idx):
            continue
        ens = row[ens_idx].strip()
        sym = row[sym_idx].strip()
        if ens and sym:
            ensembl_to_symbol[ens] = sym

# GTEx genes
gtex_path = root / 'data/gtex/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz'
gtex_genes = set()
with gzip.open(gtex_path, 'rt', encoding='utf-8', errors='ignore') as f:
    for _ in range(2):
        next(f)
    for line in f:
        parts = line.rstrip('\n').split('\t')
        if not parts:
            continue
        base = parts[0].split('.')[0]
        sym = ensembl_to_symbol.get(base)
        if sym:
            gtex_genes.add(sym)

intersection_genes = go_genes & gtex_genes

records = []
for hpo, genes in hpo_to_genes.items():
    total = len(genes)
    go_cov = len(genes & go_genes)
    gtex_cov = len(genes & gtex_genes)
    inter_cov = len(genes & intersection_genes)
    cov_frac = inter_cov/total if total else 0
    records.append((hpo, hpo_names[hpo], total, go_cov, gtex_cov, inter_cov, cov_frac))

records_20_200 = [r for r in records if 20 <= r[2] <= 200]

print(f"HP:0001083 stats: {[r for r in records if r[0]=='HP:0001083'][0] if any(r[0]=='HP:0001083' for r in records) else 'not found'}")
counts = [r[2] for r in records]
print(f"Gene count per term: min={min(counts)}, max={max(counts)}, median={median(counts)}, mean={mean(counts):.1f}")

# Largest within <=200
largest = sorted(records_20_200, key=lambda r: r[2], reverse=True)[:15]
print("\nLargest (within 20-200) top 15:")
for r in largest:
    print(f"{r[0]}\t{r[2]}\t{r[6]:.2f}\t{r[1]}")

# Mid-size around 40-120, pick best coverage then closeness to 70
mid = [r for r in records_20_200 if 40 <= r[2] <= 120]
mid_ranked = sorted(mid, key=lambda r: (r[6], -abs(r[2]-70), r[2]), reverse=True)[:20]
print("\nMid-size (40-120) top 20 by coverage then closeness to 70 genes:")
for r in mid_ranked:
    print(f"{r[0]}\t{r[2]}\t{r[6]:.2f}\t{r[1]}")

# Smaller 20-40
small = [r for r in records_20_200 if 20 <= r[2] < 40]
small_ranked = sorted(small, key=lambda r: (r[6], -r[2]), reverse=True)[:15]
print("\nSmall (20-39) top 15:")
for r in small_ranked:
    print(f"{r[0]}\t{r[2]}\t{r[6]:.2f}\t{r[1]}")
