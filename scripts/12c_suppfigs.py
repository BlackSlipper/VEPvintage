#!/usr/bin/env python
import sys; sys.path.insert(0,__import__("os").path.dirname(__import__("os").path.abspath(__file__)))
exec(open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)),"12_figures.py")).read().split("def fig1():")[0])
def figS1():
    R=pd.read_csv(f"{P}/results/dms_before_after_pub.tsv",sep="\t").dropna(subset=["auc_DMS"])
    R=R[R.assay.isin(["BRCA1_HUMAN_Findlay_2018","MSH2_HUMAN_Jia_2020"])]
    fig,axs=plt.subplots(1,2,figsize=(5.0,2.3),sharey=True)
    for ax,(a,g) in zip(axs,R.groupby("assay")):
        cols=["auc_DMS","auc_REVEL","auc_AlphaMissense","auc_ESM1b"]; names=["DMS assay","REVEL","AlphaMissense","ESM1b"]
        x=np.arange(len(cols))
        for j,(st,col) in enumerate([("before_pub",FILL["sky"]),("after_pub",FILL["rose"])]):
            r=g[g.set==st]
            if len(r)==0: continue
            ax.bar(x+(j-0.5)*0.36,r[cols].values[0],width=0.34,color=col,edgecolor=EDGE,lw=0.4,label=("Labelled up to publication year" if st=="before_pub" else "Labelled after publication")+f" (n={int(r.n.iloc[0])})")
        ax.set_xticks(x); ax.set_xticklabels(names,rotation=30,ha="right"); ax.set_ylim(0.7,1.0); ax.set_title(a.replace("_HUMAN_"," ").replace("_"," "),loc="left")
        ax.yaxis.grid(True,color=GRID,lw=0.5); ax.set_axisbelow(True); ax.legend(fontsize=5.2,loc="upper center",bbox_to_anchor=(0.5,-0.48),ncol=1)
    axs[0].set_ylabel("AUROC")
    fig.tight_layout(); save(fig,"FigS1_DMS_publication")
def figS2():
    L=pd.read_csv(f"{P}/results/leakage_by_vintage.tsv",sep="\t")
    fig,axs=plt.subplots(1,2,figsize=(6.0,2.3))
    a=axs[0]; a.bar(L.vyear,L.VARITY_insample_gain*1000,color=FILL["sky"],width=0.7,edgecolor=EDGE,lw=0.4)
    a.set_xlabel("Label vintage"); a.set_ylabel("VARITY_R minus VARITY_R_LOO\n(AUROC × 1000)"); a.set_title("In-sample advantage of VARITY",loc="left")
    a.set_xticks([2015,2018,2021,2024]); a.set_xticklabels(["≤2015","2018","2021","2024"]); a.yaxis.grid(True,color=GRID,lw=0.5); a.set_axisbelow(True)
    b=axs[1]
    for t,c in [("ClinPred",SKY),("PHACTboost",PEACH),("ESM1b",MINT)]:
        b.plot(L.vyear,L[t],marker="o",ms=2.5,lw=1.3,color=c,label=t)
    b.set_xlabel("Label vintage"); b.set_ylabel("AUROC (all variants)"); b.set_title("ClinVar-trained vs unsupervised",loc="left")
    b.set_xticks([2015,2018,2021,2024]); b.set_xticklabels(["≤2015","2018","2021","2024"]); b.legend(fontsize=5.6); b.yaxis.grid(True,color=GRID,lw=0.5); b.set_axisbelow(True)
    fig.tight_layout(); save(fig,"FigS3_leakage")
if __name__=="__main__":
    for w in sys.argv[1:]: globals()[w](); print("done",w)
