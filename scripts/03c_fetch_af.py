#!/usr/bin/env python
"""Fetch gnomAD v2.1.1 (exome/genome) AF + ExAC AF via MyVariant for ids. usage: ids out.tsv.gz shard nshard"""
import sys, gzip, time, os, requests
ids_f,out_f=sys.argv[1],sys.argv[2]; si,ns=int(sys.argv[3]),int(sys.argv[4])
ids=[l.strip() for i,l in enumerate(open(ids_f)) if l.strip() and i%ns==si]
F="gnomad_exome.af.af,gnomad_exome.af.af_popmax,gnomad_genome.af.af,gnomad_genome.af.af_popmax,dbnsfp.exac.af,gnomad_exome.an.an"
def g(d,p):
    for k in p.split("."):
        if isinstance(d,list): d=d[0] if d else None
        if not isinstance(d,dict) or k not in d: return ""
        d=d[k]
    if isinstance(d,list): d=d[0] if d else ""
    return "" if d is None else str(d)
S=requests.Session()
with gzip.open(out_f,"wt") as out:
    out.write("query\tgnomad_ex_af\tgnomad_ex_popmax\tgnomad_ge_af\tgnomad_ge_popmax\texac_af\tgnomad_ex_an\n")
    for b in range(0,len(ids),1000):
        batch=ids[b:b+1000]
        for a in range(8):
            try:
                r=S.post("https://myvariant.info/v1/variant",data={"ids":",".join(batch),"fields":F,"assembly":"hg38"},timeout=300); r.raise_for_status(); res=r.json(); break
            except Exception as e: print("retry",e,flush=True); time.sleep(10*(a+1))
        else: continue
        seen=set()
        for h in res:
            q=h.get("query")
            if q in seen: continue
            seen.add(q)
            out.write("\t".join([q]+[g(h,p) for p in ["gnomad_exome.af.af","gnomad_exome.af.af_popmax","gnomad_genome.af.af","gnomad_genome.af.af_popmax","dbnsfp.exac.af","gnomad_exome.an.an"]])+"\n")
        out.flush(); print(b+len(batch),"/",len(ids),flush=True); time.sleep(0.5)
print("DONE",flush=True)
