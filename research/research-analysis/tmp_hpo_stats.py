import gzip, csv
from collections import defaultdict
from statistics import mean
from pathlib import Path

root = Path('.')
# HPO -> genes
hpo_file = root / 'data/hpo2gene/phenotype_to_genes.txt'
hpo_to_genes = defaultdict(set)
all_hpo_genes = set()
with hpo_file.open('r', encoding='utf-8') as f:
    next(f)
    for line in f:
        hpo_id, hpo_name, ncbi, symbol, disease = line.rstrip('\n').split('\t')
        hpo_to_genes[hpo_id].add(symbol)
        all_hpo_genes.add(symbol)

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
        gene_id = parts[0]
        base = gene_id.split('.')[0]
        sym = ensembl_to_symbol.get(base)
        if sym:
            gtex_genes.add(sym)

intersection_genes = all_hpo_genes & go_genes & gtex_genes

stats = []
for hpo, genes in hpo_to_genes.items():
    total = len(genes)
    if total == 0:
        continue
    go_cov = len(genes & go_genes)
    gtex_cov = len(genes & gtex_genes)
    inter_cov = len(genes & intersection_genes)
    stats.append((hpo, total, go_cov, gtex_cov, inter_cov))

candidates = [s for s in stats if 20 <= s[1] <= 200]
candidates_ranked = sorted(
    candidates,
    key=lambda x: ((x[4]/x[1]) if x[1] else 0, -(abs(x[1]-80)), x[1]),
    reverse=True,
)

print(f"Total HPO terms: {len(stats)}")
print(f"Candidates (20-200 genes): {len(candidates)}")
print(f"Genes in any HPO term: {len(all_hpo_genes):,}")
print(f"GO genes: {len(go_genes):,}, GTEx mapped genes: {len(gtex_genes):,}, intersection gene_set-like: {len(intersection_genes):,}")

coverages = [c[4]/c[1] for c in candidates if c[1]>0]
coverages_sorted = sorted(coverages)
median_cov = coverages_sorted[len(coverages_sorted)//2] if coverages_sorted else float('nan')
print(f"Coverage (intersection) median={median_cov:.3f}, mean={mean(coverages):.3f}")

print("Top 25 candidates by coverage (hpo_id,total,inter_cov,go_cov,gtex_cov,coverage):")
for hpo, total, go_cov, gtex_cov, inter_cov in candidates_ranked[:25]:
    cov = inter_cov/total if total else 0
    print(f"{hpo}\t{total}\t{inter_cov}\t{go_cov}\t{gtex_cov}\t{cov:.2f}")

low_cov = [s for s in candidates if (s[4]/s[1]) < 0.2]
print(f"Low-coverage candidates (<0.2): {len(low_cov)}")
print("Examples (first 10):")
for hpo, total, go_cov, gtex_cov, inter_cov in low_cov[:10]:
    cov = inter_cov/total if total else 0
    print(f"{hpo}\t{total}\t{inter_cov}\t{cov:.2f}")
