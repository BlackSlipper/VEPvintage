#!/usr/bin/env python
import sys, os; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
FIG=f"{P}/figures"; os.makedirs(FIG,exist_ok=True)
plt.rcParams.update({"font.family":"Liberation Sans","font.size":7,"axes.titlesize":7.5,"axes.labelsize":7,"xtick.labelsize":6.5,"ytick.labelsize":6.5,
  "legend.fontsize":6.5,"axes.spines.top":False,"axes.spines.right":False,"axes.linewidth":0.6,"xtick.major.width":0.6,"ytick.major.width":0.6,
  "axes.edgecolor":"#444444","xtick.color":"#444444","ytick.color":"#444444","axes.labelcolor":"#1a1a1a","lines.markeredgecolor":"#4a4a4a","lines.markeredgewidth":0.4,"hatch.linewidth":0.45,"pdf.fonttype":42,"savefig.dpi":300,
  "axes.grid":False,"legend.frameon":False})
# pastel palette built around #7bd3ea, #a1eebd, #f6f7c4 and #f6d6d6
SKY="#5bb8dc"; MUSTARD="#d6bd4f"; ROSE="#e98a9e"; PERI="#8aa3e6"; MINT="#5ec88f"; LAV="#a98fd6"; PEACH="#f4a46a"; ORCHID="#d58ac9"
CAT=[SKY,MUSTARD,ROSE,PERI,MINT,LAV,PEACH,ORCHID]
FILL={"sky":"#7bd3ea","mint":"#a1eebd","butter":"#f1eca6","blush":"#f6d6d6","rose":"#f2b6c1","palesky":"#cdeef8","peach":"#fbd3b0"}
PERC=[SKY,MUSTARD,ROSE]; EDGE="#4a4a4a"; NULLC="#de6f88"
PLP_C=ROSE; BLB_C=SKY; VUS_C=MUSTARD; CONF_C=LAV
CATS=["Supervised-ClinVar","Supervised-curated","Meta-score","Proxy-label","Unsupervised-PLM","Unsupervised-MSA","Unsupervised-genomic","Conservation"]
CC=dict(zip(CATS,CAT))
GRID="#ececec"; INK2="#5c5c5c"
SEQ=LinearSegmentedColormap.from_list("pastel",["#fbfbe9","#f6f7c4","#a1eebd","#7bd3ea","#4ea6d3","#2c6f9e"])
def repel(ys,mind):
    o=np.argsort(ys); y=np.array(ys,float)[o].copy()
    for i in range(1,len(y)):
        if y[i]-y[i-1]<mind: y[i]=y[i-1]+mind
    out=np.empty_like(y); out[o]=y; return out
def lab(ax,s,y=1.06): ax.text(-0.12,y,s,transform=ax.transAxes,fontsize=9,fontweight="bold",va="bottom",ha="left")
def save(fig,name):
    fig.savefig(f"{FIG}/{name}.pdf",bbox_inches="tight"); fig.savefig(f"{FIG}/{name}.png",bbox_inches="tight",dpi=300); plt.close(fig)
META2=META.set_index("tool")
SUBA="Ambry Genetics"; SUBB="Labcorp Genetics (formerly Invitae), Labcorp"
ORDER=["<=2015Q1"]+[str(y) for y in range(2015,2027)]
XL=["≤2015Q1"]+[str(y) for y in range(2015,2027)]
REP=["REVEL","AlphaMissense","ESM1b","CADD","phyloP100way","ClinPred"]
REPC=dict(zip(REP,CAT[:6]))

def fig1():
    T=pd.read_csv(f"{P}/results/snapshot_class_counts.tsv",sep="\t")
    T["date"]=pd.to_datetime(T.Snapshot.replace("current","2026-11"))
    Tr=pd.read_csv(f"{P}/results/benchmark_tiers.tsv.gz",sep="\t")
    A=pd.read_csv(f"{P}/results/variant_level_annotated.tsv.gz",sep="\t",usecols=["hgvs38","conc_submitters"])
    Tr=Tr.merge(A.drop_duplicates("hgvs38"),on="hgvs38",how="left")
    Tr["vbin"]=np.where(Tr.vintage_first=="2015-02","<=2015Q1",Tr.vintage_first.str[:4])
    S=pd.read_csv(f"{P}/results/scv_textmining.tsv.gz",sep="\t")
    cs=S.ClinSig.str.lower(); S=S[cs.str.contains("pathogenic|benign")&~cs.str.contains("conflict|uncertain")]
    fig,ax=plt.subplots(1,3,figsize=(7.1,2.2),gridspec_kw={"width_ratios":[1,1.25,1]})
    a=ax[0]
    for c,col,nm in [("VUS",VUS_C,"VUS"),("BLB",BLB_C,"B/LB"),("PLP",PLP_C,"P/LP"),("CONF",CONF_C,"Conflicting")]:
        a.plot(T.date,T[c],color=col,lw=1.5,label=nm)
    ends=np.log10([T[c].iloc[-1] for c in ["VUS","BLB","PLP","CONF"]]); ly=repel(ends,0.3)
    for c,nm,yy in zip(["VUS","BLB","PLP","CONF"],["VUS","B/LB","P/LP","Conflicting"],ly): a.text(T.date.iloc[-1]+pd.Timedelta(days=60),10**yy,nm,va="center",fontsize=6,color="#0b0b0b")
    a.set_yscale("log"); a.set_ylabel("Missense SNVs in ClinVar"); a.set_xlim(T.date.iloc[0],T.date.iloc[-1]+pd.Timedelta(days=620)); a.set_xticks(pd.to_datetime(["2015","2018","2021","2024","2027"])); a.set_xticklabels(["2015","2018","2021","2024","2027"])
    a.set_title("ClinVar missense classifications",loc="left"); lab(a,"a")
    a.yaxis.grid(True,color=GRID,lw=0.5); a.set_axisbelow(True)
    b=ax[1]
    x=np.arange(len(ORDER))
    plp=Tr[Tr.final=="PLP"].groupby("vbin").size().reindex(ORDER).fillna(0)
    blbA=Tr[(Tr.final=="BLB")&(Tr.conc_submitters==SUBA)].groupby("vbin").size().reindex(ORDER).fillna(0)
    blbO=Tr[(Tr.final=="BLB")].groupby("vbin").size().reindex(ORDER).fillna(0)-blbA
    b.bar(x,plp/1e3,color=FILL["rose"],width=0.8,label="P/LP",edgecolor=EDGE,lw=0.4)
    b.bar(x,blbO/1e3,bottom=plp/1e3,color=FILL["sky"],width=0.8,label="B/LB (other sources)",edgecolor=EDGE,lw=0.4)
    b.bar(x,blbA/1e3,bottom=(plp+blbO)/1e3,color=FILL["palesky"],width=0.8,label="B/LB (Ambry only)",edgecolor=EDGE,lw=0.4,hatch="////")
    b.set_xticks(x); b.set_xticklabels(XL,rotation=60,ha="right"); b.set_ylabel("Variants (×1,000)"); b.set_xlabel("Label vintage (first year with final class)")
    b.set_ylim(0,40); b.legend(loc="upper left",fontsize=5.8); b.set_title("Evaluation set (n = %s) by label vintage"%format(len(Tr),","),loc="left"); lab(b,"b")
    b.yaxis.grid(True,color=GRID,lw=0.5); b.set_axisbelow(True)
    c=ax[2]
    S["year"]=pd.to_numeric(S.year,errors="coerce"); S=S[(S.year>=2013)&(S.year<=2026)]
    S["code"]=S[["PP3","BP4"]].max(axis=1)
    S["tool_named"]=S["any_tool"]
    S["generic_nonA"]=(S.any_insilico==1)&(S.Submitter!=SUBA)
    g=S.groupby("year")[["generic_nonA","tool_named","code","REVEL"]].mean()
    for col,cc,nm in [("generic_nonA",SKY,"Any computational evidence*"),("tool_named",PEACH,"Named predictor"),("code",MINT,"Explicit PP3/BP4"),("REVEL",LAV,"REVEL named")]:
        c.plot(g.index,g[col]*100,color=cc,lw=1.5,marker="o",ms=2.5,label=nm)
    c.set_ylabel("% of P/LP and B/LB SCVs"); c.set_xlabel("Date last evaluated (year)"); c.legend(loc="upper center",bbox_to_anchor=(0.5,-0.32),ncol=2,fontsize=5.8)
    c.set_title("Computational evidence cited in SCV text",loc="left"); lab(c,"c"); c.yaxis.grid(True,color=GRID,lw=0.5); c.set_axisbelow(True)
    fig.tight_layout(w_pad=1.2); save(fig,"Fig1_overview")

def fig2():
    R=pd.read_csv(f"{P}/results/auc_by_vintage_subsets.tsv",sep="\t")
    R=R[R.subset=="rare_common"]
    W=R.pivot(index="tool",columns="vbin",values="auc")[ORDER]
    cat=META2.loc[W.index,"category"]
    W=W.assign(_c=cat.map({c:i for i,c in enumerate(CATS)}),_m=W.mean(axis=1)).sort_values(["_c","_m"],ascending=[True,False])
    cats=W._c.map(dict(enumerate(CATS))); W=W.drop(columns=["_c","_m"])
    fig=plt.figure(figsize=(7.1,6.6)); G0=fig.add_gridspec(2,1,height_ratios=[2.2,1],hspace=0.62)
    top=G0[0].subgridspec(1,2,wspace=0.45); bot=G0[1].subgridspec(1,5,wspace=0.62)
    a=fig.add_subplot(top[0])
    im=a.imshow(W.values,aspect="auto",cmap=SEQ,vmin=0.65,vmax=1.0)
    a.set_yticks(range(len(W))); a.set_yticklabels(W.index,fontsize=5.4); a.set_xticks(range(len(ORDER))); a.set_xticklabels(XL,rotation=60,ha="right",fontsize=5.8)
    for i,cn in enumerate(cats): a.add_patch(plt.Rectangle((-1.15,i-0.5),0.5,1,color=CC[cn],clip_on=False,lw=0))
    a.set_xlim(-0.5,len(ORDER)-0.5); a.tick_params(axis="y",pad=14)
    cax=a.inset_axes([0.0,-0.2,0.45,0.022]); cb=fig.colorbar(im,cax=cax,orientation="horizontal"); cb.set_label("AUROC (rare variants)",fontsize=6); cb.outline.set_linewidth(0.4); cb.ax.tick_params(labelsize=5.5)
    a.set_title("AUROC by label vintage, 38 core predictors",loc="left"); lab(a,"a")
    a.spines[:].set_visible(False); a.tick_params(length=0)
    b=fig.add_subplot(top[1])
    x=np.arange(len(ORDER))
    rep=["ClinPred","REVEL","BayesDel_noAF","AlphaMissense","ESM1b","SIFT","phyloP100way"]
    ly=repel([W.loc[t].values[-1] for t in rep],0.0075)
    for t,yy in zip(rep,ly):
        y=W.loc[t].values; b.plot(x,y,color=CC[META2.loc[t,"category"]],lw=1.4,marker="o",ms=2.3)
        b.text(x[-1]+0.35,yy,t,fontsize=5.8,va="center",color="#0b0b0b")
    b.set_xticks(x); b.set_xticklabels(XL,rotation=60,ha="right",fontsize=5.8); b.set_ylabel("AUROC (rare variants)"); b.set_xlim(-0.5,len(ORDER)+2.8)
    b.yaxis.grid(True,color=GRID,lw=0.5); b.set_axisbelow(True); b.axvspan(8.5,12.5,color="#fbfbe3",zorder=0,lw=0)
    b.text(10.5,b.get_ylim()[0]+0.005,"2023–2026",ha="center",fontsize=6,color=INK2)
    b.set_title("Representative predictors",loc="left"); lab(b,"b")
    hs=[plt.Rectangle((0,0),1,1,color=CC[c]) for c in CATS]
    fig.legend(hs,CATS,loc="upper center",ncol=2,bbox_to_anchor=(0.76,0.415),fontsize=5.8,handlelength=1)
    Rb=pd.read_csv(f"{P}/results/robustness.tsv",sep="\t"); Wg=pd.read_csv(f"{P}/results/review_within_gene.tsv",sep="\t")
    per=["<=2019","2020-2022","2023-2026"]; perl=["\u22642019","2020\u201322","2023\u201326"]
    for k,t in enumerate(["REVEL","ESM1b","GERP++"]):
        c=fig.add_subplot(bot[k]); col=CC[META2.loc[t,"category"]]
        pooled=Rb[(Rb.analysis=="rare_vintage_first")&(Rb.tool==t)].set_index("period").auc.reindex(per)
        within=Wg[(Wg.analysis=="within_gene_vintage")&(Wg.tool==t)].set_index("period").auc.reindex(per)
        c.plot(range(3),pooled.values,color=col,lw=1.5,marker="o",ms=3,label="Pooled")
        c.plot(range(3),within.values,color=col,lw=1.5,marker="s",ms=3,ls="--",label="Within gene")
        c.set_xticks(range(3)); c.set_xticklabels(perl,fontsize=5.6); c.set_title(t,loc="left",fontsize=7)
        lo=min(pooled.min(),within.min()); hi=max(pooled.max(),within.max()); pad=(hi-lo)*0.35+0.005
        c.set_ylim(lo-pad,hi+pad); c.yaxis.grid(True,color=GRID,lw=0.5); c.set_axisbelow(True); c.tick_params(axis="y",labelsize=5.6)
        if k==0: c.set_ylabel("AUROC (rare)"); lab(c,"c",y=1.12); c.legend(fontsize=5.4,loc="lower left",handlelength=2.2)
    # (d) benchmarks built from successive ClinVar releases (results/r2/H_download_date_*.tsv, script 25)
    H=pd.read_csv(f"{P}/results/r2/H_download_date_auc.tsv",sep="\t"); Hn=pd.read_csv(f"{P}/results/r2/H_download_date_sets.tsv",sep="\t")
    rel=[r for r in Hn.release if Hn.set_index("release").loc[r,"n_rare_blb"]>=500]
    g=H[H.release.isin(rel)].groupby("release"); xr=np.arange(len(rel))
    d=fig.add_subplot(bot[3:5]); ink="#3d8db8"
    for col,ls,mk,name in [("pooled","-","o","Pooled"),("within","--","s","Within gene")]:
        q=g[col].quantile([0.25,0.5,0.75]).unstack().reindex(rel)
        d.plot(xr,q[0.5],color=ink,lw=1.5,ls=ls,marker=mk,ms=2.6,label=name)
    nb=Hn.set_index("release").loc[rel,"n_rare_blb"]
    i23=rel.index("2023-01"); i24=rel.index("2024-01")
    d.axvline(i23+0.5,color=INK2,lw=0.6,ls=":")
    d.text(i23+0.6,d.get_ylim()[0]+0.004,"rare B/LB %.1fk \u2192 %.1fk"%(nb.iloc[i23]/1e3,nb.iloc[i24]/1e3),fontsize=5.4,color=INK2,va="bottom")
    d.set_xticks(xr); d.set_xticklabels([("Nov 2026" if r=="2026-11" else r[:4]) for r in rel],rotation=60,ha="right",fontsize=5.6)
    d.yaxis.grid(True,color=GRID,lw=0.5); d.set_axisbelow(True); d.tick_params(axis="y",labelsize=5.6)
    d.set_ylabel("Median AUROC (rare)"); d.set_title("Benchmark built from each ClinVar release",loc="left",fontsize=7)
    d.legend(fontsize=5.4,loc="center left",bbox_to_anchor=(0.0,0.5),handlelength=2.2); lab(d,"d",y=1.12)
    save(fig,"Fig2_vintage_auc")
if __name__=="__main__":
    which=sys.argv[1:] or ["fig1","fig2"]
    for w in which: globals()[w](); print("done",w)
