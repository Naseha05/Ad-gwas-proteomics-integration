"""Map Alzheimer's GWAS lead variants (Bellenguez et al. 2022, GRCh38) to candidate genes.

Evidence layers per locus
  1. reported   : gene label printed in the paper's locus tables
  2. nearest    : nearest protein-coding gene body to the lead variant (Ensembl 111)
  3. window     : all protein-coding genes overlapping lead +/- WINDOW bp
  4. prioritized: genes with published functional prioritization (paper tier 1/2 named in text,
                  or long-standing literature consensus for known loci) -- NOT computed here.
Usage: python 01_map_loci.py <hg38.rda> <lead_variants.csv> <outdir>
"""
import sys, pandas as pd, pyreadr

WINDOW = 250_000
MERGE = 500_000  # lead variants closer than this on one chromosome form one locus

rda, leads_csv, out = sys.argv[1:4]
genes = pyreadr.read_r(rda)["hg38"]
genes = genes[genes.biotype == "protein_coding"].copy()
genes["chrom"] = genes.chrom.astype(str)
genes = genes.drop_duplicates("gene_symbol")

leads = pd.read_csv(leads_csv)
leads["chr"] = leads.chr.astype(str)
leads = leads.sort_values(["chr", "pos_grch38"]).reset_index(drop=True)

# --- merge lead variants into loci
locus_id, last = [], None
n = 0
for _, r in leads.iterrows():
    if last is None or r.chr != last.chr or r.pos_grch38 - last.pos_grch38 > MERGE:
        n += 1
    locus_id.append(f"L{n:02d}")
    last = r
leads["locus"] = locus_id

# --- literature / paper prioritization (source noted)
PAPER_TIER = {  # named as tier 1 / tier 2 in Bellenguez 2022 text
    "ICA1": 1, "ICA1L": 1, "DGKQ": 1, "DOC2A": 1, "WDR81": 1, "LIME1": 1, "EGFR": 1, "RITA1": 1,
    "MYO15A": 1, "APP": 1, "RHOH": 1, "BLNK": 1, "SIGLEC11": 1, "OTULIN": 1,
    "LILRB2": 2, "RASGEF1C": 2, "SHARPIN": 1, "RBCK1": 1,
}
# Long-standing consensus causal-gene candidates at known loci (literature, not re-derived here)
LIT = {
    "CR1": ["CR1"], "BIN1": ["BIN1"], "INPP5D": ["INPP5D"], "CLNK": ["CLNK", "HS3ST1"],
    "HLA-DQA1": ["HLA-DQA1", "HLA-DRB1"], "UNC5CL": ["TREM2", "TREML2"], "TREM2": ["TREM2"],
    "TREML2": ["TREML2"], "CD2AP": ["CD2AP"], "EPDR1": ["EPDR1"],
    "SPDYE3": ["ZCWPW1", "NYAP1", "PILRA", "PILRB"], "EPHA1": ["EPHA1"], "PTK2B": ["PTK2B"],
    "CLU": ["CLU"], "USP6NL": ["USP6NL", "ECHDC3"], "SPI1": ["SPI1", "CELF1"],
    "MS4A4A": ["MS4A4A", "MS4A6A"], "EED": ["PICALM"], "SORL1": ["SORL1"], "FERMT2": ["FERMT2"],
    "SLC24A4": ["SLC24A4", "RIN3"], "SPPL2A": ["SPPL2A"], "MINDY2": ["ADAM10"], "APH1B": ["APH1B"],
    "BCKDK": ["KAT8", "BCKDK"], "IL34": ["IL34"], "PLCG2": ["PLCG2"], "SCIMP": ["SCIMP"],
    "WNT3": ["MAPT", "WNT3"], "ABI3": ["ABI3"], "TSPOAP1": ["TSPOAP1"], "ACE": ["ACE"],
    "ABCA7": ["ABCA7"], "CASS4": ["CASS4"], "ADAMTS1": ["ADAMTS1"], "APOE": ["APOE"],
    # new loci: reported gene is also a commonly cited candidate
    "SORT1": ["SORT1"], "ADAM17": ["ADAM17"], "MME": ["MME"], "IDUA": ["IDUA"], "ANKH": ["ANKH"],
    "TNIP1": ["TNIP1"], "TMEM106B": ["TMEM106B"], "CTSB": ["CTSB"], "ABCA1": ["ABCA1"],
    "ANK3": ["ANK3"], "TSPAN14": ["TSPAN14"], "TPCN1": ["TPCN1"], "SNX1": ["SNX1"],
    "CTSH": ["CTSH"], "GRN": ["GRN"], "NCK2": ["NCK2"], "PRKD3": ["PRKD3"], "COX7C": ["COX7C"],
    "PLEKHA1": ["PLEKHA1"], "SLC2A4RG": ["SLC2A4RG"], "KLF16": ["KLF16"], "MAF": ["MAF"],
    "FOXF1": ["FOXF1"], "JAZF1": ["JAZF1"], "SEC61G": ["SEC61G", "EGFR"],
}

def dist(pos, g):
    if g.gene_start <= pos <= g.gene_end:
        return 0
    return min(abs(pos - g.gene_start), abs(pos - g.gene_end))

rows, check = [], []
for loc, grp in leads.groupby("locus", sort=False):
    ch = grp.chr.iloc[0]
    lo, hi = grp.pos_grch38.min() - WINDOW, grp.pos_grch38.max() + WINDOW
    g = genes[(genes.chrom == ch) & (genes.gene_end >= lo) & (genes.gene_start <= hi)].copy()
    lead = grp.sort_values("P", key=lambda s: pd.to_numeric(s, errors="coerce")).iloc[0]
    g["dist"] = [min(dist(p, x) for p in grp.pos_grch38) for x in g.itertuples()]
    nearest = g.sort_values("dist").gene_symbol.iloc[0] if len(g) else None
    reported = sorted(set(grp.reported_gene))
    prio = set()
    for rg in reported:
        prio.update(LIT.get(rg, []))
    for sym, t in PAPER_TIER.items():
        if sym in set(g.gene_symbol):
            prio.add(sym)
    for rg in reported:  # QC: is the reported gene inside the window?
        if rg in set(genes.gene_symbol):
            gg = genes[genes.gene_symbol == rg].iloc[0]
            d = min(dist(p, gg) for p in grp.pos_grch38) if gg.chrom == ch else -1
            check.append((loc, rg, d))
    for x in g.itertuples():
        rows.append(dict(locus=loc, chr=ch, lead_rsid=lead.rsid, lead_p=lead.P,
                         reported="/".join(reported), category=grp.category.iloc[0],
                         gene=x.gene_symbol, dist_bp=x.dist, nearest=x.gene_symbol == nearest,
                         is_reported=x.gene_symbol in reported, prioritized=x.gene_symbol in prio,
                         paper_tier=PAPER_TIER.get(x.gene_symbol, "")))
    missing = prio - set(g.gene_symbol)
    for m in missing:  # prioritized gene outside window (e.g. long-range); keep it, flagged
        rows.append(dict(locus=loc, chr=ch, lead_rsid=lead.rsid, lead_p=lead.P,
                         reported="/".join(reported), category=grp.category.iloc[0], gene=m,
                         dist_bp=None, nearest=False, is_reported=m in reported,
                         prioritized=True, paper_tier=PAPER_TIER.get(m, "")))

m = pd.DataFrame(rows)
m.to_csv(f"{out}/locus_gene_map.csv", index=False)
leads.to_csv(f"{out}/lead_variants_with_locus.csv", index=False)
qc = pd.DataFrame(check, columns=["locus", "reported_gene", "dist_bp"])
qc.to_csv(f"{out}/qc_reported_gene_distance.csv", index=False)
print("lead variants:", len(leads), "| loci:", leads.locus.nunique())
print("window genes (unique):", m[m.dist_bp.notna()].gene.nunique(),
      "| prioritized (unique):", m[m.prioritized].gene.nunique())
print("QC reported gene distance (bp) >100kb or missing:")
print(qc[(qc.dist_bp > 100_000) | (qc.dist_bp < 0)].to_string())
print("reported genes not in Ensembl protein-coding set:",
      sorted({r for r in leads.reported_gene if r not in set(genes.gene_symbol)}))
