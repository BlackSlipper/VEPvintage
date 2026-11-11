#!/usr/bin/env python
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)),"12_figures.py")).read().split("def fig1():")[0])
PER=["<=2019","2020-2022","2023-2026"]; PERL=["≤2019","2020–22","2023–26"]
def dumbbell(ax,A,B,labA,labB,colA,colB,title,leg="top"):
    d=pd.DataFrame({"a":A,"b":B}).dropna().sort_values("b")
    y=np.arange(len(d))
    ax.hlines(y,d.b,d.a,color="#cfcfcf",lw=0.8,zorder=1)
    ax.scatter(d.b,y,s=10,color=colB,zorder=2,label=labB,edgecolor=EDGE,lw=0.35)
    ax.scatter(d.a,y,s=10,color=colA,zorder=3,label=labA,edgecolor=EDGE,lw=0.35)
    ax.set_yticks(y); ax.set_yticklabels(d.index,fontsize=5.4); ax.set_ylim(-0.8,len(d)-0.2)
    ax.xaxis.grid(True,color=GRID,lw=0.5); ax.set_axisbelow(True)
    if leg=="top":
        ax.legend(loc="lower center",bbox_to_anchor=(0.4,1.02),fontsize=5.6,handletextpad=0.2,ncol=1); ax.set_title(title,loc="left",pad=24)
    else:
        ax.legend(loc="upper center",bbox_to_anchor=(0.45,-0.09),fontsize=5.6,handletextpad=0.2); ax.set_title(title,loc="left")
    return d
def fig3():
    E=pd.read_csv(f"{P}/results/evidence_source_contrasts.tsv",sep="\t")
    V=pd.read_csv(f"{P}/results/variant_level_annotated.tsv.gz",sep="\t",usecols=["final","rare","vyear","legacy","conc_submitters"])
    fig=plt.figure(figsize=(7.1,6.6)); gs=fig.add_gridspec(2,3,height_ratios=[1,2.1],hspace=0.38,wspace=0.6)
    # b: AUROC by period and benign-label source (rare)
    for k,t in enumerate(["REVEL","phyloP100way"]):
        a=fig.add_subplot(gs[0,k+1])
        for j,(key,nm,col) in enumerate([("BLB_all","All B/LB",SKY),("BLB_excluding_Ambry_Labcorp","B/LB excl. Ambry, Labcorp",MINT),("BLB_only_from:Ambry Genetics","B/LB from Ambry only",MUSTARD)]):
            ys=[];xs=[]
            for i,p in enumerate(PER):
                r=E[(E.contrast==f"rare:{p}:{key}")&(E.tool==t)]
                if len(r): xs.append(i+(j-1)*0.12); ys.append(r.auc.iloc[0])
            a.plot(xs,ys,marker="o",ms=3.5,lw=1.2,color=col,label=nm)
        a.set_xticks(range(3)); a.set_xticklabels(PERL); a.set_ylabel("AUROC (rare variants)"); a.set_title(t,loc="left")
        a.yaxis.grid(True,color=GRID,lw=0.5); a.set_axisbelow(True); a.set_xlim(-0.4,2.4)
        if k==0: lab(a,"b"); lo_,hi_=a.get_ylim(); a.set_ylim(lo_,hi_+(hi_-lo_)*0.38); a.legend(loc="upper left",fontsize=5.5)
    # a: share of rare B/LB by source per vintage year
    b=fig.add_subplot(gs[0,0])
    W=V[(V.final=="BLB")&(V.rare)].copy(); W["vb"]=np.where(W.legacy,2015,W.vyear.clip(lower=2015))
    W["src"]=np.where(W.conc_submitters==SUBA,"Ambry only",np.where(W.conc_submitters==SUBB,"Labcorp only","Other / multiple"))
    ct=pd.crosstab(W.vb,W.src); ct=ct[["Other / multiple","Labcorp only","Ambry only"]]
    bot=np.zeros(len(ct))
    for c,col in zip(ct.columns,[FILL["sky"],FILL["mint"],FILL["butter"]]):
        b.bar(ct.index,ct[c]/1e3,bottom=bot,color=col,width=0.8,edgecolor=EDGE,lw=0.4,label=c); bot+=ct[c].values/1e3
    b.set_ylabel("Rare B/LB variants (×1,000)"); b.set_xlabel("Label vintage"); b.set_title("Source of rare B/LB labels",loc="left"); lab(b,"a")
    b.legend(loc="upper left",fontsize=5.5); b.yaxis.grid(True,color=GRID,lw=0.5); b.set_axisbelow(True)
    b.set_xticks([2015,2018,2021,2024]); b.set_xticklabels(["≤2015","2018","2021","2024"])
    # c,d: explicit evidence codes (rare)
    c=fig.add_subplot(gs[1,0:1]); d=fig.add_subplot(gs[1,1:2])
    piv=E[E.contrast.str.startswith("rare:")&~E.contrast.str.contains(":20|:<=")].pivot(index="tool",columns="contrast",values="auc")
    dumbbell(c,piv["rare:BLB_with_BP4"],piv["rare:BLB_codes_without_BP4"],"B/LB citing BP4","B/LB citing other codes only",ROSE,SKY,"Benign side: BP4")
    c.set_xlabel("AUROC (vs all rare P/LP)"); lab(c,"c",y=1.17)
    dumbbell(d,piv["rare:PLP_with_PP3"],piv["rare:PLP_codes_without_PP3"],"P/LP citing PP3","P/LP citing other codes only",ROSE,SKY,"Pathogenic side: PP3")
    d.set_xlabel("AUROC (vs all rare B/LB)"); lab(d,"d",y=1.17)
    e=fig.add_subplot(gs[1,2])
    dl=pd.DataFrame({"BP4":piv["rare:BLB_with_BP4"]-piv["rare:BLB_codes_without_BP4"],"PP3":piv["rare:PLP_with_PP3"]-piv["rare:PLP_codes_without_PP3"]})
    dl["cat"]=META2.loc[dl.index,"category"]
    for cn in CATS:
        s=dl[dl.cat==cn]; e.scatter(s.BP4,s.PP3,s=14,color=CC[cn],edgecolor=EDGE,lw=0.35,label=cn)
    for t in ["ESM1b","PrimateAI","REVEL","GERP++","ClinPred","AlphaMissense"]:
        if t in dl.index: e.annotate(t,(dl.loc[t,"BP4"],dl.loc[t,"PP3"]),fontsize=5.3,xytext=(3,2),textcoords="offset points")
    e.set_xlabel("AUROC inflation, BP4-cited B/LB"); e.set_ylabel("AUROC inflation, PP3-cited P/LP"); e.set_title("Inflation by predictor",loc="left"); lab(e,"e")
    e.xaxis.grid(True,color=GRID,lw=0.5); e.yaxis.grid(True,color=GRID,lw=0.5); e.set_axisbelow(True)
    e.legend(loc="upper center",bbox_to_anchor=(0.45,-0.16),ncol=2,fontsize=5.2,handletextpad=0.1,columnspacing=0.6)
    save(fig,"Fig3_provenance")
    dl.to_csv(f"{P}/results/evidence_code_inflation.tsv",sep="\t")
def fig4():
    B=pd.read_csv(f"{P}/results/benchmark_full_vs_clean.tsv",sep="\t")
    O=pd.read_csv(f"{P}/results/relative_slopes.tsv",sep="\t")
    Pr=pd.read_csv(f"{P}/results/prospective_auc.tsv",sep="\t")
    core=O[O.subset=="rare_common"].tool.tolist()
    fig=plt.figure(figsize=(7.1,4.8)); gs=fig.add_gridspec(1,3,width_ratios=[1.0,0.9,0.9],wspace=1.05)
    a=fig.add_subplot(gs[0])
    B2=pd.read_csv(f"{P}/results/benchmark_clean_rare_recent.tsv",sep="\t")
    W=B2[B2.tool.isin(core)].pivot(index="tool",columns="set",values="auc")
    rf=W.full_rare_2023_26.rank(ascending=False); rc=W.clean_rare_2023_26.rank(ascending=False)
    for t in W.index:
        col=CC[META2.loc[t,"category"]]; a.plot([0,1],[rf[t],rc[t]],color=col,lw=1.0,marker="o",ms=2.5)
        a.text(-0.04,rf[t],t,ha="right",va="center",fontsize=4.8); a.text(1.04,rc[t],t,ha="left",va="center",fontsize=4.8)
    a.set_ylim(len(W)+0.5,0.5); a.set_xlim(-0.05,1.05); a.set_xticks([0,1]); a.set_xticklabels(["All labels","Provenance-\nfiltered"]); a.set_yticks([])
    a.spines["left"].set_visible(False); a.set_title("Rank, rare variants labelled 2023–26",loc="left"); lab(a,"a")
    b=fig.add_subplot(gs[1])
    o=O[O.subset=="rare_common"].sort_values("rel_slope")
    y=np.arange(len(o))
    b.hlines(y,o.rel_slope-1.96*o.se,o.rel_slope+1.96*o.se,color="#cfcfcf",lw=0.8)
    b.scatter(o.rel_slope,y,s=12,c=[CC[c] for c in o.category],edgecolor=EDGE,lw=0.35,zorder=3)
    b.axvline(0,color=INK2,lw=0.6); b.set_yticks(y); b.set_yticklabels(o.tool,fontsize=5.0)
    b.set_xlabel("Relative vintage slope\n(Δ logit-AUROC vs median per year)"); b.set_title("Who gains over vintages",loc="left"); lab(b,"b")
    b.xaxis.grid(True,color=GRID,lw=0.5); b.set_axisbelow(True)
    c=fig.add_subplot(gs[2])
    p=Pr[Pr["T"]=="2024-01"].pivot(index="tool",columns="set",values="auc").loc[core]; p.index=[NM(t) for t in p.index]
    dumbbell(c,p["known_at_T"],p["VUS_resolved_after_T"],"Labelled before 2024","VUS resolved after 2024-01",SKY,ROSE,"Time-split benchmark",leg="inside")
    c.set_xlabel("AUROC"); lab(c,"c")
    save(fig,"Fig4_rankings")
def fig5():
    fig,axs=plt.subplots(2,2,figsize=(7.1,4.4),sharey="row")
    tools=["REVEL","BayesDel_noAF","VEST4","MutPred2","AlphaMissense","VARITY_R"]
    for col,(fn,ttl) in enumerate([("calibration_LR_by_period_rare.tsv","All rare variants"),("calibration_LR_by_period_rare_clean.tsv","Provenance-filtered rare variants")]):
        if not os.path.exists(f"{P}/results/{fn}"): continue
        R=pd.read_csv(f"{P}/results/{fn}",sep="\t")
        for row,ev in enumerate(["PP3_Supporting","PP3_Moderate"]):
            a=axs[row,col]
            for j,p in enumerate(PER):
                r=R[(R.evidence==ev)&(R.period==p)].set_index("tool").reindex(tools)
                x=np.arange(len(tools))+(j-1)*0.22
                a.errorbar(x,r.LR,yerr=[r.LR-r.lo,r.hi-r.LR],fmt="o",ms=3,color=PERC[j],elinewidth=0.7,capsize=0,label=PERL[j])
            a.axhline(R[R.evidence==ev].nominal.iloc[0],color=INK2,ls="--",lw=0.7)
            a.text(len(tools)-0.5,R[R.evidence==ev].nominal.iloc[0]*1.08,"nominal OddsPath",fontsize=5.3,ha="right",color=INK2)
            a.set_yscale("log"); from matplotlib.ticker import FixedLocator,NullFormatter,FormatStrFormatter
            tk=[2,3,4,6,8,10] if row==0 else [5,10,20,40]
            a.yaxis.set_major_locator(FixedLocator(tk)); a.yaxis.set_major_formatter(FormatStrFormatter("%g")); a.yaxis.set_minor_formatter(NullFormatter())
            a.set_xticks(range(len(tools))); a.set_xticklabels(tools if row==1 else [],rotation=35,ha="right")
            a.set_title(ev.replace(chr(95),chr(32))+" — "+ttl,loc="left"); a.yaxis.grid(True,color=GRID,lw=0.5,which="both"); a.set_axisbelow(True)
            if col==0: a.set_ylabel("Likelihood ratio (P/LP : B/LB)")
            if row==0 and col==0: a.legend(title="Label vintage",fontsize=5.6,title_fontsize=5.8,loc="upper left",ncol=3)
    lab(axs[0,0],"a"); lab(axs[0,1],"b")
    fig.tight_layout(); save(fig,"Fig5_calibration")
NM=lambda t: t.replace("100way","").replace("PolyPhen2_","PolyPhen2-")
def fig4r2():
    A2=pd.read_csv(f"{P}/results/r2/A_full_minus_filtered.tsv",sep="\t")
    RK=pd.read_csv(f"{P}/results/r2/B_ranks_core_common.tsv",sep="\t")
    Pr=pd.read_csv(f"{P}/results/prospective_auc.tsv",sep="\t")
    core=RK.tool.unique().tolist()
    fig=plt.figure(figsize=(7.1,5.0)); gs=fig.add_gridspec(1,3,width_ratios=[1.0,0.95,0.9],wspace=1.15)
    a=fig.add_subplot(gs[0])
    d=A2[(A2["filter"]=="clean")&(A2.period=="2023-2026")&A2.tool.isin(core)].set_index("tool").sort_values("d_pooled")
    y=np.arange(len(d))
    a.hlines(y-0.15,d.d_pooled_lo,d.d_pooled_hi,color=ROSE,lw=0.7,alpha=0.8); a.scatter(d.d_pooled,y-0.15,s=9,color=ROSE,zorder=3,label="Pooled",edgecolor=EDGE,lw=0.3)
    a.hlines(y+0.15,d.d_within_lo,d.d_within_hi,color=SKY,lw=0.7,alpha=0.8); a.scatter(d.d_within,y+0.15,s=9,color=SKY,zorder=3,label="Within gene",edgecolor=EDGE,lw=0.3)
    a.axvline(0,color=INK2,lw=0.6); a.set_yticks(y); a.set_yticklabels([NM(t) for t in d.index],fontsize=5.0); a.set_ylim(-0.8,len(d)-0.2)
    a.set_xlabel("AUROC, all minus filtered labels"); a.xaxis.grid(True,color=GRID,lw=0.5); a.set_axisbelow(True)
    a.legend(loc="lower right",fontsize=5.6,handletextpad=0.2,frameon=False); a.set_title("Drop on filtering, rare 2023–26",loc="left"); lab(a,"a")
    b=fig.add_subplot(gs[1])
    X=RK[RK.metric=="pooled"].pivot(index="tool",columns="set",values="rank")
    for t in X.index:
        col=CC[META2.loc[t,"category"]]; b.plot([0,1],[X.loc[t,"full"],X.loc[t,"clean"]],color=col,lw=1.0,marker="o",ms=2.5)
        b.text(-0.04,X.loc[t,"full"],NM(t),ha="right",va="center",fontsize=4.8); b.text(1.04,X.loc[t,"clean"],NM(t),ha="left",va="center",fontsize=4.8)
    b.set_ylim(len(X)+0.5,0.5); b.set_xlim(-0.05,1.05); b.set_xticks([0,1]); b.set_xticklabels(["All labels","Provenance-\nfiltered"]); b.set_yticks([])
    b.spines["left"].set_visible(False); b.set_title("Rank, 38 core predictors",loc="left"); lab(b,"b")
    c=fig.add_subplot(gs[2])
    p=Pr[Pr["T"]=="2024-01"].pivot(index="tool",columns="set",values="auc").loc[core]; p.index=[NM(t) for t in p.index]
    dumbbell(c,p["known_at_T"],p["VUS_resolved_after_T"],"Labelled before 2024","VUS resolved after 2024-01",SKY,ROSE,"Time-split benchmark",leg="inside")
    c.set_xlabel("AUROC"); lab(c,"c")
    save(fig,"Fig4_rankings")
def fig5r2():
    from matplotlib.ticker import FixedLocator,NullFormatter,FormatStrFormatter
    L=pd.read_csv(f"{P}/results/r2/C_LR_pooled_MH.tsv",sep="\t")
    tools=["REVEL","BayesDel_noAF","VEST4","MutPred2","CADD","PrimateAI","MPC","phyloP100way","AlphaMissense","ESM1b","VARITY_R"]
    lbl={"BayesDel_noAF":"BayesDel","phyloP100way":"phyloP","VARITY_R":"VARITY_R"}
    fig,axs=plt.subplots(2,3,figsize=(7.1,4.6),sharey="row")
    nom={"PP3_Supporting":2.08,"PP3_Moderate":4.33}
    for col,(st,est,ttl) in enumerate([("full","pooled","All labels, pooled"),("full","MH","All labels, within gene"),("clean","pooled","Filtered labels, pooled")]):
        for row,ev in enumerate(["PP3_Supporting","PP3_Moderate"]):
            a=axs[row,col]
            for j,p in enumerate(PER):
                r=L[(L.evidence==ev)&(L.set==st)&(L.period==p)].set_index("tool").reindex(tools)
                v,lo,hi=(r.LR_pooled,r.LRp_lo,r.LRp_hi) if est=="pooled" else (r.LR_MH,r.LRmh_lo,r.LRmh_hi)
                v=v.astype(float); lo=lo.astype(float).clip(lower=0.3); hi=hi.astype(float).clip(upper=200)
                x=np.arange(len(tools))+(j-1)*0.24
                a.errorbar(x,v,yerr=[np.maximum(v-lo,0),np.maximum(hi-v,0)],fmt="o",ms=2.4,color=PERC[j],elinewidth=0.6,capsize=0,label=PERL[j])
            a.axhline(nom[ev],color=INK2,ls="--",lw=0.6)
            a.set_yscale("log")
            tk=[1,2,4,8,16] if row==0 else [2,5,10,20,50]
            a.set_ylim((0.7,24) if row==0 else (1.5,90))
            a.yaxis.set_major_locator(FixedLocator(tk)); a.yaxis.set_major_formatter(FormatStrFormatter("%g")); a.yaxis.set_minor_formatter(NullFormatter())
            a.set_xticks(range(len(tools))); a.set_xticklabels([lbl.get(t,t) for t in tools] if row==1 else [],rotation=50,ha="right",fontsize=5.2)
            a.set_title((ev.replace(chr(95),chr(32))+"\n"+ttl),loc="left",fontsize=6.4); a.yaxis.grid(True,color=GRID,lw=0.5,which="major"); a.set_axisbelow(True)
            if col==0: a.set_ylabel("Likelihood ratio (P/LP : B/LB)")
            if row==0 and col==0: a.legend(title="Label vintage",fontsize=5.4,title_fontsize=5.6,loc="upper left",ncol=3,handletextpad=0.1,columnspacing=0.5)
    lab(axs[0,0],"a"); lab(axs[0,1],"b"); lab(axs[0,2],"c")
    fig.tight_layout(); save(fig,"Fig5_calibration")
def figS3():
    O=pd.read_csv(f"{P}/results/relative_slopes_rare_fdr.tsv",sep="\t").sort_values("rel_slope")
    fig,b=plt.subplots(figsize=(3.4,5.0))
    y=np.arange(len(O))
    b.hlines(y,O.rel_slope-1.96*O.se,O.rel_slope+1.96*O.se,color="#cfcfcf",lw=0.8)
    b.scatter(O.rel_slope,y,s=12,c=[CC[c] for c in O.category],edgecolor=EDGE,lw=0.35,zorder=3)
    b.axvline(0,color=INK2,lw=0.6); b.set_yticks(y); b.set_yticklabels(O.tool,fontsize=5.0)
    b.set_xlabel("Relative vintage slope\n(Δ logit-AUROC vs median per year)"); b.set_title("Relative change over vintages",loc="left")
    b.xaxis.grid(True,color=GRID,lw=0.5); b.set_axisbelow(True)
    fig.tight_layout(); save(fig,"FigS2_relative_slopes")
if __name__=="__main__":
    for w in sys.argv[1:]: globals()[w](); print("done",w)
