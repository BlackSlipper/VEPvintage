#!/usr/bin/env python
"""Within-gene and paired analyses on rare variants (max gnomAD v2/ExAC AF < 1e-4):
(A) pooled vs within-gene AUROC for full and provenance-filtered tiers by period; paired gene-cluster bootstrap (B=1000);
    decomposition of the full-minus-filtered difference into within-gene and between-gene parts; gene-matched variant.
(B) ranks of the 38 core predictors on the common subset (rare 2023-26), pooled and within-gene, bootstrap rank intervals.
(C) Mantel-Haenszel within-gene likelihood ratios for ClinGen score intervals by period (vs pooled LRs).
(D) diagnosis of the filtered drop on <=2019 labels: removal reasons, evaluation dates of flagging records, as-of-vintage filter.
(E) vintage jump 2023-26 minus <=2019: paired bootstrap, pooled and within-gene; equivalence margin 0.02.
(F) resolved VUS at T=2024-01: CIs for provenance-filtered subgroups."""
import sys, os; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
from scipy.stats import kendalltau, rankdata
from multiprocessing import Pool
B=int(__import__("os").environ.get("VEP_BOOT","1000")); OUT=f"{P}/results/r2"; os.makedirs(OUT,exist_ok=True)
pd.set_option("display.width",260); pd.set_option("display.max_rows",600); pd.set_option("display.max_columns",40)
D=load()
T=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t",usecols=["hgvs38","clean","rare","period","explicit_comp","comp_nonambry","ambry_only"])
D=D.merge(T,on="hgvs38")
A_=pd.read_csv(f"{P}/results/variant_level_annotated2.tsv.gz",sep="\t",usecols=["hgvs38","conc_submitters"]).drop_duplicates("hgvs38")
D=D.merge(A_,on="hgvs38",how="left")
for c in ["clean","rare","explicit_comp","comp_nonambry","ambry_only"]: D[c]=D[c].fillna(False).astype(bool)
SA="Ambry Genetics"; SB="Labcorp Genetics (formerly Invitae), Labcorp"
S=pd.read_csv(f"{P}/results/scv_textmining.tsv.gz",sep="\t",dtype={"VariationID":str},usecols=["VariationID","SCV","Submitter","ClinSig","has_text","any_insilico","year"])
C=pd.read_csv(f"{P}/results/scv_acmg_codes.tsv.gz",sep="\t",dtype={"VariationID":str})
st=S[S.has_text==1].groupby("Submitter").agg(n_text=("SCV","size"),frac_comp=("any_insilico","mean"))
cc=C.assign(comp=C.codes.str.contains(r"\bPP3\b|\bBP4\b")).groupby("Submitter").comp.mean().rename("fcc")
st=st.join(cc).fillna({"fcc":0}); noncomp=set(st[(st.n_text>=20)&(st.frac_comp<0.05)&(st.fcc<0.05)].index)-{SA}
cs=S.ClinSig.str.lower(); S["cl"]=np.where(cs.str.contains("pathogenic"),"PLP",np.where(cs.str.contains("benign"),"BLB","O"))
S=S.merge(D[["VariationID","final","vyear"]].drop_duplicates("VariationID"),on="VariationID"); S=S[S.cl==S.final]
D["subclean"]=D.VariationID.isin(set(S[S.Submitter.isin(noncomp)].VariationID))
D["Bonly"]=D.conc_submitters.fillna("")==SB
tools=[t for t in TOOLS if t in D.columns]; cov=D[tools].notna().mean(); core=[t for t in tools if cov[t]>=0.90]
print("n eval",len(D),"rare",D.rare.sum(),"core",len(core),"noncomp submitters",len(noncomp))
R=D[D.rare].copy()
genes=pd.Index(sorted(D.Gene.dropna().unique())); ng=len(genes)
rng=np.random.default_rng(2026); W=np.stack([np.bincount(rng.integers(0,ng,ng),minlength=ng) for _ in range(B)]).astype(float)
def gstats(G,t):
    g=G[["Gene","y",t]].dropna(); gi=genes.get_indexer(g.Gene.values)
    r=g.groupby("Gene")[t].rank().values; y=g.y.values.astype(float)
    n1=np.bincount(gi,weights=y,minlength=ng); n=np.bincount(gi,minlength=ng).astype(float)
    rs=np.bincount(gi,weights=r*y,minlength=ng)
    return rs-n1*(n1+1)/2, n1*(n-n1)
def pprep(G,t):
    g=G[["Gene","y",t]].dropna(); lev=pd.Series(g[t].values).rank(method="dense").values.astype(int)-1
    return lev,lev.max()+1,g.y.values.astype(float),genes.get_indexer(g.Gene.values)
def pauc(pr,w):
    lev,K,y,gi=pr; ww=w[gi]
    Pk=np.bincount(lev,weights=ww*y,minlength=K); Nk=np.bincount(lev,weights=ww*(1-y),minlength=K)
    cN=np.cumsum(Nk)-Nk; return (Pk*(cN+0.5*Nk)).sum()/(Pk.sum()*Nk.sum())
ONE=np.ones(ng)
def ci(x,q=(2.5,97.5)): x=np.asarray(x,float); x=x[np.isfinite(x)]; return np.percentile(x,q) if len(x) else (np.nan,np.nan)
per_list=["<=2019","2020-2022","2023-2026"]
SETS={}
for per in per_list:
    G=R[R.period==per]
    SETS[("full",per)]=G; SETS[("clean",per)]=G[G.clean]; SETS[("subclean",per)]=G[G.subclean]; SETS[("clean_exclBonly",per)]=G[G.clean&~G.Bonly]
def jobA(k):
    s,p,t=k; G=SETS[(s,p)]; pr=pprep(G,t)
    bt=np.array([pauc(pr,W[b]) for b in range(B)]); U,pa=gstats(G,t)
    return k,dict(n=len(pr[2]),n_plp=int(pr[2].sum()),genes_with_pairs=int((pa>0).sum()),pooled=pauc(pr,ONE),within=U.sum()/pa.sum()),bt,(W@U)/(W@pa),U,pa
KEYS=[(s,p,t) for (s,p) in SETS for t in core]
with Pool(20) as pool: res=pool.map(jobA,KEYS,chunksize=4)
BT={};WB={};UP={};PT={};rows=[]
for k,info,bt,wb,U,pa in res:
    BT[k]=bt;WB[k]=wb;UP[k]=(U,pa);PT[k]=(info["pooled"],info["within"])
    a,b=ci(bt);c,d=ci(wb);e,f=ci(bt-wb)
    rows.append(dict(set=k[0],period=k[1],tool=k[2],**info,pooled_lo=a,pooled_hi=b,within_lo=c,within_hi=d,between=info["pooled"]-info["within"],between_lo=e,between_hi=f))
A1=pd.DataFrame(rows); A1.to_csv(f"{OUT}/A_tiers_pooled_within.tsv",sep="\t",index=False)
print("=== A1 median over core: pooled / within by set x period ===")
print(A1.groupby(["set","period"])[["n","n_plp","genes_with_pairs","pooled","within","between"]].median().round(3).to_string())
key=[t for t in ["REVEL","AlphaMissense","ESM1b","CADD","GERP++","MetaRNN","BayesDel_noAF","VEST4","MutPred2","phyloP100way","PrimateAI","SIFT"] if t in core]
print("key tools in core:",key)
print(A1[A1.tool.isin(key)&(A1.period=="2023-2026")].pivot(index="tool",columns="set",values=["pooled","within"]).round(3).to_string())
rows=[]
for p in per_list:
    for x in ["clean","subclean","clean_exclBonly"]:
        for t in core:
            f=("full",p,t); c=(x,p,t)
            dp=PT[f][0]-PT[c][0]; dw=PT[f][1]-PT[c][1]; bp=BT[f]-BT[c]; bw=WB[f]-WB[c]
            Uf,pf=UP[f]; Uc,pc=UP[c]; m=pc>0
            dm=Uf[m].sum()/pf[m].sum()-PT[c][1]; bm=(W[:,m]@Uf[m])/(W[:,m]@pf[m])-WB[c]
            a,b=ci(bp);c1,d1=ci(bw);e,f1=ci(bm)
            rows.append(dict(period=p,filter=x,tool=t,d_pooled=dp,d_pooled_lo=a,d_pooled_hi=b,d_within=dw,d_within_lo=c1,d_within_hi=d1,
                             d_within_genematched=dm,d_wgm_lo=e,d_wgm_hi=f1,frac_within=dw/dp if abs(dp)>1e-9 else np.nan))
A2=pd.DataFrame(rows); A2.to_csv(f"{OUT}/A_full_minus_filtered.tsv",sep="\t",index=False)
print("=== A2 full minus filtered: medians over core; n tools with within CI>0 ===")
g=A2.groupby(["filter","period"])
print(pd.concat([g[["d_pooled","d_within","d_within_genematched"]].median(),g.apply(lambda d:(d.d_within_lo>0).sum()).rename("n_within_CI_gt0"),
                 g.apply(lambda d:(d.d_within_hi<0).sum()).rename("n_within_CI_lt0"),g.d_within.min().rename("dw_min"),g.d_within.max().rename("dw_max"),
                 (g.d_within.sum()/g.d_pooled.sum()).rename("share_within")],axis=1).round(3).to_string())
print(A2[A2.tool.isin(key)&(A2.period=="2023-2026")&(A2["filter"]=="clean")].round(3).to_string())
print(A2[A2.tool.isin(key)&(A2.period=="2023-2026")&(A2["filter"]=="subclean")].round(3).to_string())
print(A2[A2.tool.isin(key)&(A2.period=="<=2019")&(A2["filter"]=="clean")].round(3).to_string())
# (E) vintage jump
rows=[]
for x in ["full","clean","subclean"]:
    for t in core:
        a=(x,"2023-2026",t); b=(x,"<=2019",t)
        bp=BT[a]-BT[b]; bw=WB[a]-WB[b]
        rows.append(dict(set=x,tool=t,jump_pooled=PT[a][0]-PT[b][0],jp_lo=ci(bp)[0],jp_hi=ci(bp)[1],jump_within=PT[a][1]-PT[b][1],
                         jw_lo=ci(bw)[0],jw_hi=ci(bw)[1],jw_lo90=ci(bw,(5,95))[0],jw_hi90=ci(bw,(5,95))[1]))
E=pd.DataFrame(rows); E["equiv02"]=(E.jw_lo90>-0.02)&(E.jw_hi90<0.02); E.to_csv(f"{OUT}/E_jump.tsv",sep="\t",index=False)
print("=== E jump (2023-26 minus <=2019) ===")
g=E.groupby("set")
print(pd.concat([g[["jump_pooled","jump_within"]].median(),g.apply(lambda d:(d.jp_lo>0).sum()).rename("n_pooled_CI_gt0"),
      g.apply(lambda d:(d.jw_lo>0).sum()).rename("n_within_CI_gt0"),g.apply(lambda d:(d.jw_hi<0).sum()).rename("n_within_CI_lt0"),g.equiv02.sum().rename("n_equiv_0.02")],axis=1).round(3).to_string())
print(E[E.tool.isin(key)&(E.set=="full")].round(3).to_string())
# (B) ranks on common subset
Cc=R[(R.period=="2023-2026")&R[core].notna().all(axis=1)]
RS={"full":Cc,"clean":Cc[Cc.clean],"subclean":Cc[Cc.subclean]}
print("common subset rare 2023-26:",{k:(len(v),int(v.y.sum())) for k,v in RS.items()})
def jobB(k):
    s,t=k; G=RS[s]; pr=pprep(G,t); U,pa=gstats(G,t)
    return k,pauc(pr,ONE),np.array([pauc(pr,W[b]) for b in range(B)]),U.sum()/pa.sum(),(W@U)/(W@pa)
with Pool(20) as pool: rb=pool.map(jobB,[(s,t) for s in RS for t in core],chunksize=2)
rows=[]; MAT={}
for (s,t),pp,bp,ww,bw in rb:
    rows.append(dict(set=s,tool=t,pooled=pp,within=ww)); MAT.setdefault((s,"pooled"),{})[t]=bp; MAT.setdefault((s,"within"),{})[t]=bw
Rk=pd.DataFrame(rows)
out=[]
for s in RS:
    for m in ["pooled","within"]:
        M=np.column_stack([MAT[(s,m)][t] for t in core]); RkB=rankdata(-M,axis=1)
        pt=Rk[Rk.set==s].set_index("tool").loc[core,m]
        out.append(pd.DataFrame({"set":s,"metric":m,"tool":core,"auc":pt.values,"rank":rankdata(-pt.values),
                                 "rank_lo":np.percentile(RkB,2.5,axis=0),"rank_hi":np.percentile(RkB,97.5,axis=0)}))
RK=pd.concat(out); RK.to_csv(f"{OUT}/B_ranks_core_common.tsv",sep="\t",index=False)
print("=== B ranks (38 core, common subset, rare 2023-26) ===")
for m in ["pooled","within"]:
    X=RK[RK.metric==m]; Wr=X.pivot(index="tool",columns="set",values="rank"); Lo=X.pivot(index="tool",columns="set",values="rank_lo"); Hi=X.pivot(index="tool",columns="set",values="rank_hi")
    Wr["d_clean"]=Wr.full-Wr.clean; Wr["d_sub"]=Wr.full-Wr.subclean
    Wr["robust_up"]=(Hi.clean<Lo.full)&(Hi.subclean<Lo.full); Wr["robust_down"]=(Lo.clean>Hi.full)&(Lo.subclean>Hi.full)
    Wr["robust_up_clean_only"]=(Hi.clean<Lo.full); Wr["robust_down_clean_only"]=(Lo.clean>Hi.full)
    for c in ["full","clean","subclean"]: Wr[c+"_CI"]=Lo[c].round().astype(int).astype(str)+"-"+Hi[c].round().astype(int).astype(str)
    print("--",m); print(Wr.sort_values("clean").to_string())
    A=X.pivot(index="tool",columns="set",values="auc")
    print(m,"Kendall tau full-clean %.3f full-subclean %.3f clean-subclean %.3f"%(kendalltau(A.full,A.clean)[0],kendalltau(A.full,A.subclean)[0],kendalltau(A.clean,A.subclean)[0]),
          "| max |d| clean",Wr.d_clean.abs().max(),"sub",Wr.d_sub.abs().max(),"| robust up",list(Wr.index[Wr.robust_up]),"down",list(Wr.index[Wr.robust_down]))
PW=RK.pivot_table(index="tool",columns=["set","metric"],values="auc")
print("tau pooled vs within: full %.3f clean %.3f subclean %.3f"%tuple(kendalltau(PW[(s,"pooled")],PW[(s,"within")])[0] for s in ["full","clean","subclean"]))
# (C) Mantel-Haenszel within-gene LRs
Th=pd.read_csv(f"{P}/meta/clingen_thresholds.tsv",sep="\t"); Th["lo"]=Th.lo.astype(float); Th["hi"]=Th.hi.astype(float)
Th=Th[Th.rawcol.isin(D.columns)]; print("calibrated tools with scores:",sorted(Th.tool.unique()),len(Th.tool.unique()))
def lrstats(G,rc,lo,hi):
    g=G[["Gene","y",rc]].dropna(); y=g.y.values.astype(float); inb=((g[rc]>=lo)&(g[rc]<hi)).values.astype(float); gi=genes.get_indexer(g.Gene.values)
    a=np.bincount(gi,weights=inb*y,minlength=ng); b=np.bincount(gi,weights=inb*(1-y),minlength=ng)
    n1=np.bincount(gi,weights=y,minlength=ng); n0=np.bincount(gi,weights=1-y,minlength=ng); N=np.maximum(n1+n0,1); ok=(n1>0)&(n0>0)
    return a,b,n1,n0,np.where(ok,a*n0/N,0),np.where(ok,b*n1/N,0)
def pooledLR(a,b,n1,n0): return ((a+.5)/(n1+.5))/((b+.5)/(n0+.5))
LSETS={}
for x in ["full","clean","subclean"]:
    for per in per_list: LSETS[(x,per)]=SETS[(x,per)]
LSETS[("subclean","<=2022")]=R[R.subclean&(R.period!="2023-2026")]; LSETS[("clean","<=2022")]=R[R.clean&(R.period!="2023-2026")]; LSETS[("full","<=2022")]=R[R.period!="2023-2026"]
rows=[]; LB={}
for _,r in Th.iterrows():
    inv=r.evidence.startswith("BP4")
    for (x,per),G in LSETS.items():
        a,b,n1,n0,num,den=lrstats(G,r.rawcol,r.lo,r.hi)
        pl=pooledLR(a.sum(),b.sum(),n1.sum(),n0.sum()); plb=pooledLR(W@a,W@b,W@n1,W@n0)
        with np.errstate(divide="ignore",invalid="ignore"):
            mh=num.sum()/den.sum() if den.sum()>0 else np.nan; mhb=(W@num)/(W@den)
        if inv: pl,plb,mh,mhb=1/pl,1/plb,1/mh,1/mhb
        LB[(r.tool,r.evidence,x,per)]=(plb,mhb)
        rows.append(dict(tool=r.tool,evidence=r.evidence,set=x,period=per,n_plp=int(n1.sum()),n_blb=int(n0.sum()),in_plp=int(a.sum()),in_blb=int(b.sum()),
                         genes_both=int(((n1>0)&(n0>0)).sum()),LR_pooled=pl,LRp_lo=ci(plb)[0],LRp_hi=ci(plb)[1],LR_MH=mh,LRmh_lo=ci(mhb)[0],LRmh_hi=ci(mhb)[1]))
L=pd.DataFrame(rows)
rows=[]
for (tool,ev),_ in L.groupby(["tool","evidence"]):
    for x,base in [("full","<=2019"),("clean","<=2019"),("subclean","<=2022"),("subclean","<=2019")]:
        if (tool,ev,x,"2023-2026") not in LB: continue
        p1,m1=LB[(tool,ev,x,"2023-2026")]; p0,m0=LB[(tool,ev,x,base)]
        q=L[(L.tool==tool)&(L.evidence==ev)&(L.set==x)].set_index("period")
        rp=q.loc["2023-2026","LR_pooled"]/q.loc[base,"LR_pooled"]; rm=q.loc["2023-2026","LR_MH"]/q.loc[base,"LR_MH"]
        rows.append(dict(tool=tool,evidence=ev,set=x,base=base,ratio_pooled=rp,rp_lo=ci(p1/p0)[0],rp_hi=ci(p1/p0)[1],ratio_MH=rm,rmh_lo=ci(m1/m0)[0],rmh_hi=ci(m1/m0)[1]))
LR=pd.DataFrame(rows); L.to_csv(f"{OUT}/C_LR_pooled_MH.tsv",sep="\t",index=False); LR.to_csv(f"{OUT}/C_LR_ratios.tsv",sep="\t",index=False)
print("=== C PP3_Supporting / BP4_Supporting LRs ===")
for ev in ["PP3_Supporting","BP4_Supporting","PP3_Moderate"]:
    q=L[L.evidence==ev]; print("--",ev); print(q.pivot_table(index="tool",columns=["set","period"],values=["LR_pooled","LR_MH"]).round(2).to_string())
    print(LR[LR.evidence==ev].round(2).to_string())
# (D) <=2019 diagnosis
Cd=C.copy(); Cd["yr"]=pd.to_datetime(Cd.Date,errors="coerce").dt.year
Cd=Cd[((Cd.final=="PLP")&Cd.codes.str.contains(r"\bPP3\b"))|((Cd.final=="BLB")&Cd.codes.str.contains(r"\bBP4\b"))]
St=S[(S.Submitter!=SA)&(S.any_insilico==1)][["VariationID","year"]].rename(columns={"year":"yr"})
F=pd.concat([Cd[["VariationID","yr"]].assign(src="code"),St.assign(src="text")])
F=F.merge(R[["VariationID","vyear","period"]].drop_duplicates("VariationID"),on="VariationID")
F["early"]=F.yr.isna()|(F.yr<=F.vyear)
R["flag_any"]=R.explicit_comp|R.comp_nonambry
R["flag_early"]=R.VariationID.map(F.groupby("VariationID").early.max()).fillna(False).astype(bool)
R["flag_minyr"]=R.VariationID.map(F.groupby("VariationID").yr.min())
R["clean_asof"]=(~R.ambry_only)&(~R.flag_early)
print("=== D consistency: flag_any vs in F ===",pd.crosstab(R.flag_any,R.VariationID.isin(set(F.VariationID))).to_string())
print("flagging-record years by period:\n",F.groupby("period").yr.describe().round(1).to_string())
print("share of flagging records evaluated after vintage year, by period:",F.groupby("period").apply(lambda d:(~d.early).mean()).round(3).to_dict())
rows=[]
for per in ["<=2019","2023-2026"]:
    G=R[R.period==per]
    for cls in ["PLP","BLB"]:
        g=G[G.final==cls]
        rows.append(dict(period=per,cls=cls,n=len(g),legacy=g.legacy.mean(),code=g.explicit_comp.mean(),text=g.comp_nonambry.mean(),A_only=g.ambry_only.mean(),
                         removed=(~g.clean).mean(),flag_any=g.flag_any.mean(),flag_late_only=(g.flag_any&~g.flag_early).mean(),removed_asof=(~g.clean_asof).mean()))
    sub={"full":G,"clean":G[G.clean],"clean_asof":G[G.clean_asof],"no_code":G[~G.explicit_comp],"no_text":G[~G.comp_nonambry],"no_Aonly":G[~G.ambry_only],
         "no_code_no_text":G[~G.explicit_comp&~G.comp_nonambry],"legacy_only":G[G.legacy],"nonlegacy":G[~G.legacy],"clean_nonlegacy":G[G.clean&~G.legacy]}
    rr=[]
    for nm,g in sub.items():
        for t in key:
            U,pa=gstats(g,t); x=g[[t,"y"]].dropna()
            rr.append(dict(period=per,subset=nm,tool=t,n=len(x),n_plp=int(x.y.sum()),pooled=auc(x.y.values,x[t].values.astype(float)),within=U.sum()/pa.sum() if pa.sum()>0 else np.nan))
    Q=pd.DataFrame(rr); Q.to_csv(f"{OUT}/D_subsets_{per.replace('<=','le')}.tsv",sep="\t",index=False)
    print("--",per); print(Q.pivot(index="subset",columns="tool",values="pooled").round(3).to_string()); print(Q.pivot(index="subset",columns="tool",values="within").round(3).to_string())
    print(Q.groupby("subset")[["n","n_plp"]].max().to_string())
Dd=pd.DataFrame(rows); Dd.to_csv(f"{OUT}/D_removal_reasons.tsv",sep="\t",index=False); print(Dd.round(3).to_string())
# (F) resolved VUS at T=2024-01
def state_at(h,Tt):
    st_="ABSENT"
    for tok in h.split(";"):
        s,c_=tok.split(":")
        if s=="current" or s>Tt: break
        st_=c_
    return st_
D["stT"]=D.history.fillna("").map(lambda h: state_at(h,"2024-01") if h else "ABSENT")
V=D[D.stT.isin(["VUS","CONF"])]
rows=[]
for nm,G in {"resolved_all":V,"resolved_clean":V[V.clean],"resolved_nonclean":V[~V.clean],"resolved_subclean":V[V.subclean]}.items():
    for t in key:
        pr=pprep(G,t); bt=[pauc(pr,W[b]) for b in range(B)]; lo,hi=ci(bt)
        rows.append(dict(set=nm,tool=t,n=len(pr[2]),n_plp=int(pr[2].sum()),auc=pauc(pr,ONE),lo=lo,hi=hi))
Fr=pd.DataFrame(rows); Fr.to_csv(f"{OUT}/F_resolved_VUS.tsv",sep="\t",index=False); print("=== F ===\n",Fr.round(3).to_string())
print("DONE")
