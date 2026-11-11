#!/usr/bin/env python
"""Genes behind the composition effect: rare B/LB labels of 2023-2026 in genes without any rare P/LP label in the
evaluation set; concentration across genes and sources (Submitter A only, Submitter B only, other)."""
import sys, os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from common import *
OUT=f"{P}/results/r2"; pd.set_option("display.width",220)
T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","Gene","final","period","rare","clean"])
A=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",usecols=["hgvs38","conc_submitters"]).drop_duplicates("hgvs38")
D=T.merge(A,on="hgvs38",how="left"); D=D[D.rare.fillna(False).astype(bool)]
SA="Ambry Genetics"; SB="Labcorp Genetics (formerly Invitae), Labcorp"
D["src"]=np.where(D.conc_submitters.fillna("")==SA,"A_only",np.where(D.conc_submitters.fillna("")==SB,"B_only","other"))
plp_genes=set(D[D.final=="PLP"].Gene.dropna())
b=D[(D.final=="BLB")&(D.period=="2023-2026")].copy(); b["no_plp_gene"]=~b.Gene.isin(plp_genes)
x=b[b.no_plp_gene]
print("rare B/LB 2023-26:",len(b),"in genes without rare P/LP:",len(x),"(%.1f%%)"%(100*len(x)/len(b)),"genes:",x.Gene.nunique(),"of",b.Gene.nunique())
g=x.groupby("Gene").agg(n_blb=("hgvs38","size"),A_only=("src",lambda s:(s=="A_only").mean()),B_only=("src",lambda s:(s=="B_only").mean()),
                       clean=("clean",lambda s:s.fillna(False).astype(bool).mean())).sort_values("n_blb",ascending=False)
cum=g.n_blb.cumsum()/g.n_blb.sum()
print("genes needed for 50%% / 80%% of these labels: %d / %d"%((cum<0.5).sum()+1,(cum<0.8).sum()+1))
print("median labels per gene %.0f; genes with >=10 labels %d"%(g.n_blb.median(),(g.n_blb>=10).sum()))
print("source of these labels:",x.src.value_counts(normalize=True).round(3).to_dict(),"| kept by per-variant filter: %.3f"%x.clean.fillna(False).astype(bool).mean())
print(g.head(25).round(2).to_string())
y=b[~b.no_plp_gene]
print("for comparison, genes with rare P/LP: source",y.src.value_counts(normalize=True).round(3).to_dict())
g.reset_index().to_csv(f"{OUT}/I_gene_drivers_blb_2023_26.tsv",sep="\t",index=False)
