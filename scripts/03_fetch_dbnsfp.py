#!/usr/bin/env python
"""Fetch dbNSFP 4.8a annotations via MyVariant.info batch API.
usage: 03_fetch_dbnsfp.py ids.txt out.tsv.gz [shard_index n_shards]
ids are hg38 HGVS genomic ids like chr17:g.7675088C>T
"""
import sys, json, gzip, time, os, requests
ids_f, out_f = sys.argv[1], sys.argv[2]
si, ns = (int(sys.argv[3]), int(sys.argv[4])) if len(sys.argv)>4 else (0,1)
ids=[l.strip() for l in open(ids_f) if l.strip()]
ids=[x for i,x in enumerate(ids) if i%ns==si]
RANK = {
 "SIFT":"sift.converted_rankscore","SIFT4G":"sift4g.converted_rankscore",
 "PolyPhen2_HDIV":"polyphen2.hdiv.rankscore","PolyPhen2_HVAR":"polyphen2.hvar.rankscore",
 "LRT":"lrt.converted_rankscore","MutationTaster":"mutationtaster.converted_rankscore",
 "MutationAssessor":"mutationassessor.rankscore","FATHMM":"fathmm.converted_rankscore",
 "PROVEAN":"provean.converted_rankscore","VEST4":"vest4.rankscore","MetaSVM":"metasvm.rankscore",
 "MetaLR":"metalr.rankscore","MetaRNN":"metarnn.rankscore","M-CAP":"m-cap.rankscore",
 "REVEL":"revel.rankscore","MutPred2":"mutpred.rankscore","MVP":"mvp.rankscore","gMVP":"gmvp.rankscore",
 "MPC":"mpc.rankscore","PrimateAI":"primateai.rankscore","DEOGEN2":"deogen2.rankscore",
 "BayesDel_addAF":"bayesdel.add_af.rankscore","BayesDel_noAF":"bayesdel.no_af.rankscore",
 "ClinPred":"clinpred.rankscore","LIST-S2":"list-s2.rankscore",
 "VARITY_R":"varity.r.rankscore","VARITY_ER":"varity.er.rankscore",
 "VARITY_R_LOO":"varity.r_loo.rankscore","VARITY_ER_LOO":"varity.er_loo.rankscore",
 "ESM1b":"esm1b.rankscore","EVE":"eve.rankscore","AlphaMissense":"alphamissense.rankscore",
 "MutFormer":"mutformer.rankscore","PHACTboost":"phactboost.rankscore",
 "CADD":"cadd.raw_rankscore","DANN":"dann.rankscore",
 "FATHMM-MKL":"fathmm-mkl.coding_rankscore","FATHMM-XF":"fathmm-xf.coding_rankscore",
 "Eigen":"eigen.raw_coding_rankscore","Eigen-PC":"eigen-pc.raw_coding_rankscore",
 "GERP++":"gerp++.rs_rankscore","phyloP100way":"phylop.100way_vertebrate.rankscore",
 "phastCons100way":"phastcons.100way_vertebrate.rankscore","SiPhy":"siphy_29way.logodds_rankscore",
}
RAW = {  # raw scores; multi-transcript -> most damaging
 "REVEL_raw":("revel.score","max"),"BayesDel_noAF_raw":("bayesdel.no_af.score","max"),
 "AlphaMissense_raw":("alphamissense.score","max"),"CADD_phred":("cadd.phred","max"),
 "VEST4_raw":("vest4.score","max"),"MutPred2_raw":("mutpred.score","max"),
 "ESM1b_raw":("esm1b.score","min"),"PrimateAI_raw":("primateai.score","max"),
 "MPC_raw":("mpc.score","max"),"PolyPhen2_HVAR_raw":("polyphen2.hvar.score","max"),
 "SIFT_raw":("sift.score","min"),"phyloP100way_raw":("phylop.100way_vertebrate.score","max"),
 "VARITY_R_raw":("varity.r.score","max"),"EVE_raw":("eve.score","max"),"MetaRNN_raw":("metarnn.score","max"),
 "ClinPred_raw":("clinpred.score","max"),
}
def get(d,path):
    for k in path.split("."):
        if isinstance(d,dict) and k in d: d=d[k]
        else: return None
    return d
def agg(v,how):
    if v is None: return None
    if isinstance(v,list):
        v=[x for x in v if isinstance(x,(int,float))]
        if not v: return None
        return max(v) if how=="max" else min(v)
    return v if isinstance(v,(int,float)) else None
def rk(v):
    if isinstance(v,list):
        v=[x for x in v if isinstance(x,(int,float))]
        return max(v) if v else None
    return v if isinstance(v,(int,float)) else None
cols=["query","found","uniprot","aapos","aaref","aaalt","genename"]+list(RANK)+list(RAW)
done=set()
if os.path.exists(out_f):
    with gzip.open(out_f,"rt") as fh:
        next(fh,None)
        for l in fh: done.add(l.split("\t",1)[0])
todo=[x for x in ids if x not in done]
print(f"shard {si}/{ns}: {len(ids)} ids, {len(done)} done, {len(todo)} todo", flush=True)
fields="dbnsfp"
S=requests.Session()
mode="at" if done else "wt"
with gzip.open(out_f,mode) as out:
    if not done: out.write("\t".join(cols)+"\n")
    for b in range(0,len(todo),1000):
        batch=todo[b:b+1000]
        for attempt in range(8):
            try:
                r=S.post("https://myvariant.info/v1/variant",data={"ids":",".join(batch),"fields":fields,"assembly":"hg38"},timeout=300)
                r.raise_for_status(); res=r.json(); break
            except Exception as e:
                print("retry",attempt,e,flush=True); time.sleep(10*(attempt+1))
        else:
            print("FAILED batch",b,flush=True); continue
        seen=set()
        for h in res:
            q=h.get("query")
            if q in seen: continue
            seen.add(q)
            n=h.get("dbnsfp") if not h.get("notfound") else None
            if n is None:
                out.write("\t".join([q,"0"]+[""]*(len(cols)-2))+"\n"); continue
            if isinstance(n,list): n=n[0]
            if not isinstance(n,dict):
                out.write("\t".join([q,"0"]+[""]*(len(cols)-2))+"\n"); continue
            up=get(n,"uniprot"); acc=""
            if isinstance(up,list): acc=",".join(sorted({u.get("acc","") for u in up if isinstance(u,dict)}))
            elif isinstance(up,dict): acc=up.get("acc","") if isinstance(up.get("acc"),str) else ",".join(up.get("acc") or [])
            aa=n.get("aa",{}) or {}
            if isinstance(aa,list):
                pl=[];rf="";al=""
                for a in aa:
                    if isinstance(a,dict):
                        p=a.get("pos"); pl+= p if isinstance(p,list) else [p]
                        rf=rf or str(a.get("ref","")); al=al or str(a.get("alt",""))
                aa={"pos":[p for p in pl if p is not None],"ref":rf,"alt":al}
            pos=aa.get("pos"); pos=",".join(map(str,sorted(set(pos)))) if isinstance(pos,list) else str(pos or "")
            gn=n.get("genename"); gn=",".join(sorted(set(gn))) if isinstance(gn,list) else str(gn or "")
            row=[q,"1",acc,pos,str(aa.get("ref","")),str(aa.get("alt","")),gn]
            for k,p in RANK.items():
                v=rk(get(n,p)); row.append("" if v is None else f"{v:.6g}")
            for k,(p,how) in RAW.items():
                v=agg(get(n,p),how); row.append("" if v is None else f"{v:.6g}")
            out.write("\t".join(row)+"\n")
        out.flush()
        print(f"{b+len(batch)}/{len(todo)}",flush=True)
        time.sleep(1)
print("DONE",flush=True)
