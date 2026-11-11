#!/usr/bin/env python
"""gnomAD v4.1 joint (exome+genome) AC/AN/FAF95-max for evaluation variants, queried per gene from the gnomAD GraphQL API."""
import pandas as pd, requests, time, os
P=__import__("os").environ.get("VEP_ROOT",__import__("os").path.dirname(__import__("os").path.dirname(__import__("os").path.abspath(__file__))))
T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","Gene"])
T["vid"]=T.hgvs38.str.replace("chr","",regex=False).str.replace(":g.","-",regex=False).str.replace(r"(\d+)([ACGT])>([ACGT])$",r"\1-\2-\3",regex=True)
want=set(T.vid); genes=sorted(set(g for gs in T.Gene.dropna() for g in str(gs).split(";")))
out=f"{P}/results/gnomad4_af.tsv"; dfile=f"{P}/results/gnomad4_genes_done.txt"
done=set(open(dfile).read().split()) if os.path.exists(dfile) else set()
q="query($g:String!){ gene(gene_symbol:$g, reference_genome:GRCh38){ variants(dataset:gnomad_r4){ variant_id joint{ ac an fafmax{ faf95_max } } } } }"
S=requests.Session(); S.headers["User-Agent"]="academic-research-script"
fo=open(out,"a"); fd=open(dfile,"a")
print(len(genes),"genes;",len(done),"done",flush=True)
for i,g in enumerate(genes):
    if g in done: continue
    d=None
    for a in range(6):
        try:
            r=S.post("https://gnomad.broadinstitute.org/api",json={"query":q,"variables":{"g":g}},timeout=180)
            if r.status_code==429:
                time.sleep(30*(a+1)); continue
            d=r.json(); break
        except Exception as e:
            time.sleep(10*(a+1))
    if d and d.get("data") and d["data"].get("gene"):
        for v in d["data"]["gene"]["variants"]:
            j=v.get("joint")
            if v["variant_id"] in want and j:
                f=(j.get("fafmax") or {}).get("faf95_max")
                fo.write("\t".join([v["variant_id"],str(j.get("ac")),str(j.get("an")),"" if f is None else str(f)])+"\n")
    fd.write(g+"\n"); fo.flush(); fd.flush()
    if i%200==0: print(i,g,flush=True)
    time.sleep(0.4)
print("DONE",flush=True)
