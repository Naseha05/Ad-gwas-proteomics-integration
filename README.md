# Alzheimer's GWAS x brain proteomics — integrated analysis

Run from this folder (Python 3, pandas, scipy, statsmodels, pyreadr):

    python scripts/01_map_loci.py data/ensembl111_hg38_genes.rda data/bellenguez2022_lead_variants.csv results
    python scripts/02_integrate.py . data/ensembl111_hg38_genes.rda

## data/
- bellenguez2022_lead_variants.csv — 83 lead variants (Bellenguez et al. 2022, Tables 1-2, GRCh38) + APOE rs429358
- ensembl111_hg38_genes.rda — Ensembl 111 gene coordinates (CRAN package enshuman)
- hgnc_symbol_ensembl_babelgene.csv — symbol -> Ensembl ID (CRAN package babelgene)
- agora_raw_lookup.txt — per-gene AMP-AD Agora proteomics results (TMT/LFQ/SRM, AD vs control), read 2026-10-17
- agora_lookup_list.csv — genes queried (89 GWAS-prioritized + 40 random background, seed 2026)
- *.gmt — GO BP 2021 (Enrichr), KEGG (MSigDB 7.5.1), Hallmark (MSigDB 7.0)

## results/
- locus_gene_map.csv — every locus x gene with distance, nearest/reported/prioritized flags
- qc_reported_gene_distance.csv — position sanity check
- proteomics_lookup_parsed.csv — parsed proteomics per gene
- overlap_tests.csv — detection and DA rate vs background (Fisher)
- integrated_gene_table.csv — gene-level evidence tiers (A-D) with proteomics and pathways
- pathway_enrichment_gwas_prioritized.csv / pathway_enrichment_window_genes.csv

To replace the per-gene lookups with the full proteome, download "TMT Proteomics for Agora" from Synapse
and swap it in for agora_raw_lookup.txt; the background test then uses all quantified proteins.
