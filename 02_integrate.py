"""Integrate GWAS candidate genes with AD brain proteomics (AMP-AD Agora lookups) and pathways.
Usage: python 02_integrate.py <project_dir>
"""
import re, sys, itertools, numpy as np, pandas as pd
from scipy.stats import fisher_exact, hypergeom, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import pyreadr

D = sys.argv[1]
ALPHA = 0.05

# ---------- parse proteomics lookups
rows = []
bg_section = False
for line in open(f"{D}/data/agora_raw_lookup.txt"):
    line = line.strip()
    if not line:
        continue
    if line.startswith("#"):
        if "random background" in line:
            bg_section = True
        continue
    sym = line.split("|")[0].strip()
    rec = dict(gene=sym, set="random_background" if bg_section else "gwas_prioritized",
               in_agora="NOT_FOUND" not in line)
    parts = [p.strip() for p in line.split("|")[1:]]
    # TMT (take the first / canonical entry)
    tmt = next((p for p in parts if p.startswith("TMT")), "TMT=NA")
    m = re.search(r"logfc=([-\deE.]+) adjp=([-\deE.]+)", tmt)
    rec["tmt_log2fc"] = float(m.group(1)) if m else np.nan
    rec["tmt_fdr"] = float(m.group(2)) if m else np.nan
    # LFQ / SRM: record minimum FDR entry and its sign
    for plat in ("LFQ", "SRM"):
        seg = next((p for p in parts if p.startswith(plat)), f"{plat}=NA")
        hits = re.findall(r"([A-Za-z]+[^ ]*) logfc=([-\deE.]+) adjp=([-\deE.]+)", seg)
        if hits:
            best = min(hits, key=lambda h: float(h[2]))
            rec[f"{plat.lower()}_best_region"] = best[0]
            rec[f"{plat.lower()}_best_log2fc"] = float(best[1])
            rec[f"{plat.lower()}_best_fdr"] = float(best[2])
            rec[f"{plat.lower()}_n_sig"] = sum(float(h[2]) < ALPHA for h in hits)
        else:
            rec[f"{plat.lower()}_best_fdr"] = np.nan
            rec[f"{plat.lower()}_n_sig"] = 0
    rows.append(rec)
P = pd.DataFrame(rows)
P["detected_tmt"] = P.tmt_fdr.notna()
P["detected_any"] = P[["tmt_fdr", "lfq_best_fdr", "srm_best_fdr"]].notna().any(axis=1)
P["da_tmt"] = P.tmt_fdr < ALPHA
P["da_any"] = (P.tmt_fdr < ALPHA) | (P.lfq_n_sig > 0) | (P.srm_best_fdr < ALPHA)
P.to_csv(f"{D}/results/proteomics_lookup_parsed.csv", index=False)

# ---------- GWAS gene-level annotation
M = pd.read_csv(f"{D}/results/locus_gene_map.csv")
pri = M[M.prioritized].drop_duplicates("gene")[["gene", "locus", "lead_rsid", "lead_p", "reported",
                                                 "category", "dist_bp", "paper_tier"]]
G = pri.merge(P[P.set == "gwas_prioritized"], on="gene", how="left")

# ---------- 1. detection & DA rate: GWAS prioritized vs random background
out = []
def rate_test(col, sub=None):
    a = P[P.set == "gwas_prioritized"]; b = P[(P.set == "random_background") & P.in_agora]
    if sub:
        a, b = a[a[sub]], b[b[sub]]
    t = [[a[col].sum(), (~a[col]).sum()], [b[col].sum(), (~b[col]).sum()]]
    orr, p = fisher_exact(t)
    out.append(dict(test=f"{col}" + (f" | {sub}" if sub else ""), gwas_k=t[0][0], gwas_n=len(a),
                    gwas_frac=t[0][0] / max(len(a), 1), bg_k=t[1][0], bg_n=len(b),
                    bg_frac=t[1][0] / max(len(b), 1), fisher_OR=orr, fisher_p=p))
rate_test("detected_tmt"); rate_test("detected_any")
rate_test("da_tmt", "detected_tmt"); rate_test("da_any", "detected_any")
R = pd.DataFrame(out)
R.to_csv(f"{D}/results/overlap_tests.csv", index=False)

# ---------- 2. pathway enrichment of GWAS prioritized genes (protein-coding background)
genes = pyreadr.read_r(sys.argv[2])["hg38"] if len(sys.argv) > 2 else None
universe = set(genes[genes.biotype == "protein_coding"].gene_symbol) if genes is not None else None

def read_gmt(path, src):
    gs = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        gs[f"{src}: {f[0]}"] = {g.split(",")[0] for g in f[2:] if g}
    return gs
GS = {}
GS.update(read_gmt(f"{D}/data/enrichr.GO_Biological_Process_2021.gmt", "GO_BP"))
GS.update(read_gmt(f"{D}/data/c2.cp.kegg.v7.5.1.symbols.gmt", "KEGG"))
GS.update(read_gmt(f"{D}/data/h.all.v7.0.symbols.gmt", "Hallmark"))
U = universe & set().union(*GS.values())
query = set(pri.gene) & U
DA = set(P[(P.set == "gwas_prioritized") & P.da_any].gene)
er = []
for name, s in GS.items():
    s = s & U
    if not 5 <= len(s) <= 500:
        continue
    k = len(query & s)
    if k < 3:
        continue
    p = hypergeom.sf(k - 1, len(U), len(s), len(query))
    er.append(dict(pathway=name, set_size=len(s), overlap=k, expected=len(query) * len(s) / len(U),
                   p=p, genes=",".join(sorted(query & s)),
                   da_protein_genes=",".join(sorted(query & s & DA))))
E = pd.DataFrame(er)
E["fdr"] = multipletests(E.p, method="fdr_bh")[1]
E["fold"] = E.overlap / E.expected
E = E.sort_values("p")
E.to_csv(f"{D}/results/pathway_enrichment_gwas_prioritized.csv", index=False)

# ---------- 3. robustness: pathway enrichment using the broader +/-250kb window genes
W = set(M[M.dist_bp.notna()].gene) & U
ew = []
for name, s in GS.items():
    s = s & U
    if not 5 <= len(s) <= 500:
        continue
    k = len(W & s)
    if k < 3:
        continue
    ew.append(dict(pathway=name, overlap=k, p=hypergeom.sf(k - 1, len(U), len(s), len(W)),
                   loci=len({l for l, g in M[M.dist_bp.notna()][["locus", "gene"]].values if g in s})))
EW = pd.DataFrame(ew); EW["fdr"] = multipletests(EW.p, method="fdr_bh")[1]
EW.sort_values("p").to_csv(f"{D}/results/pathway_enrichment_window_genes.csv", index=False)

# ---------- 4. integrated evidence table
def tier(r):
    if r.da_any and r.da_tmt:
        return "A: DA in TMT (DLPFC, ~400 brains) + GWAS"
    if r.da_any:
        return "B: DA in LFQ/SRM only + GWAS"
    if r.detected_any:
        return "C: detected, not DA"
    if r.in_agora:
        return "D: not quantified in brain proteomics"
    return "D: not quantified in brain proteomics"
G["evidence_tier"] = G.apply(tier, axis=1)
gene2paths = {}
top = E[E.fdr < 0.05]
for _, r in top.iterrows():
    for g in r.genes.split(","):
        gene2paths.setdefault(g, []).append(r.pathway)
G["top_enriched_pathways"] = G.gene.map(lambda g: "; ".join(gene2paths.get(g, [])[:3]))
G = G.sort_values(["evidence_tier", "tmt_fdr"])
G.to_csv(f"{D}/results/integrated_gene_table.csv", index=False)

print(R.to_string(index=False))
print("\nEvidence tiers:\n", G.evidence_tier.value_counts().to_string())
print("\nTier A/B genes:\n", G[G.evidence_tier.str[0].isin(["A", "B"])][
    ["gene", "locus", "tmt_log2fc", "tmt_fdr", "lfq_best_region", "lfq_best_log2fc", "lfq_best_fdr",
     "srm_best_log2fc", "srm_best_fdr"]].to_string(index=False))
print("\nTop pathways (prioritized genes): universe", len(U), "query", len(query))
print(E.head(25)[["pathway", "overlap", "fold", "p", "fdr", "genes", "da_protein_genes"]].to_string(index=False))
print("\nTop pathways (window genes): n", len(W))
print(EW.sort_values("p").head(10).to_string(index=False))
