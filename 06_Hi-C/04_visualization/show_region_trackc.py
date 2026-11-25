#!/usr/bin/env python3
# env: hic-viz-env
# REQUIRES (set in config.sh):
#   HI_C_DIR           output root of the Hi-C module (loop/, TAD/, AB_compartments/, figures/)
#   CHROMHMM_SEGMENTS  chromHMM 10-state segments BED
#   HISTONE_BW_DIR     bigwig directory <mark>/<mark>_FE.bw
#   RNA_BW_CRO         RNA-seq bigwig, control (Cro)
#   GENOME_GFF         gene annotation GFF3
#   GENOME_BED12       gene annotation in BED12 format
#   EXPR_DETAIL        optional expression table
#   HIC_FIG_DIR        optional, default <HI_C_DIR>/figures
import os, re, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pyBigWig

HI_C_DIR = os.environ.get("HI_C_DIR", os.path.join(os.getcwd(), "hi_c"))
OUTDIR = os.environ.get("HIC_FIG_DIR", os.path.join(HI_C_DIR, "figures"))
os.makedirs(OUTDIR, exist_ok=True)
SEG   = os.environ.get("CHROMHMM_SEGMENTS", "")
BWD   = os.environ.get("HISTONE_BW_DIR", "")
RNASEQ= os.environ.get("RNA_BW_CRO", "")
GFF   = os.environ.get("GENOME_GFF", "")
EXPR  = os.environ.get("EXPR_DETAIL", "")   # optional
BED12 = os.environ.get("GENOME_BED12", "")

# (tag, chrom, start, end, expand_kb, fill_color | None, use_trackc_gene)
REGIONS = [
    ("scaffold4_6363458-6663458_APP",     "scaffold_4", 6_363_458, 6_663_458, 300, "#1F77B4", True),
    ("scaffold4_193160225-193460225_ANXA2","scaffold_4",193_160_225,193_460_225, 0, None, False),
    ("scaffold2_299718914-300018914_EZR", "scaffold_2",299_718_914,300_018_914, 0, None, False),
]

state_names = {1:"ReprPC",2:"Quies",3:"Heterochr",4:"Tx",5:"Tx+Enh",
               6:"EnhW",7:"EnhA",8:"Prom",9:"PromW",10:"PromEnh"}
state_colors = {1:"#1F78B4",2:"#A6CEE3",3:"#CAB2D6",4:"#33A02C",5:"#B2DF8A",
                6:"#EF99BF",7:"#E31A1C",8:"#FF7F00",9:"#FFFF99",10:"#E9C60D"}
marks_colors = {"H3K4me3":"#E41A1C","H3K27ac":"#FF7F00","H3K4me1":"#FFD200",
                "H3K36me3":"#4DAF4A","H3K9me3":"#555555","H3K27me3":"#984EA3"}
MARK_ORDER = ["H3K4me3","H3K27ac","H3K4me1","H3K36me3","H3K9me3","H3K27me3"]

def load_rna_expr():
    expr={}
    if not EXPR or not os.path.exists(EXPR):
        return expr
    for i,l in enumerate(open(EXPR)):
        if i==0: continue
        f=l.rstrip().split("\t")
        expr[f[0]]={"exp":f[5],"maxfpkm":f[4],"meanC":f[10] if f[10]!="NA" else "NA",
                    "change":f[9] if f[9]!="NA" else "NA"}
    return expr
expr=load_rna_expr()

def load_segments(chr_,s0,e0):
    segs=[]
    for l in open(SEG):
        p=l.split()
        if p[0]!=chr_: continue
        s,e,st=int(p[1]),int(p[2]),int(p[3][1:])
        if e>s0 and s<e0: segs.append((max(s,s0),min(e,e0),st))
    segs.sort(); return segs

def load_genes(chr_,s0,e0):
    tx_gene={}; tx_strand={}; tx_name={}; tx_exons={}
    for l in open(GFF):
        if l.startswith("#"): continue
        p=l.rstrip().split("\t")
        if len(p)<9: continue
        chrom,ftype,s,e,strand,attrs=p[0],p[2],int(p[3]),int(p[4]),p[6],p[8]
        if chrom!=chr_ or e<s0 or s>e0: continue
        if ftype in ("mRNA","lnc_RNA"):
            m=re.search(r"ID=([^;]+)",attrs); g=re.search(r"geneID=([^;]+)",attrs)
            mid=m.group(1) if m else None
            tx_gene[mid]=g.group(1) if g else mid
            tx_strand[mid]=strand; tx_exons.setdefault(mid,[])
            if g:
                pn=re.search(r"EvidenceProteinID=[^;]*:([^;]+?)\s*\[",attrs)
                if pn: tx_name.setdefault(g.group(1),pn.group(1))
        elif ftype=="exon":
            m=re.search(r"Parent=([^;]+)",attrs)
            mid=m.group(1) if m else None
            if mid in tx_exons: tx_exons[mid].append((s,e))
    gene_exons={}; gene_strand={}
    for mid,exs in tx_exons.items():
        if not exs: continue
        gid=tx_gene.get(mid,mid)
        gene_exons.setdefault(gid,[]).extend(exs)
        gene_strand[gid]=tx_strand.get(mid,"+")
    merged={}
    for gid,exs in gene_exons.items():
        exs=sorted(set(exs)); m=[]
        for a,b in exs:
            if m and a<=m[-1][1]+50: m[-1]=(m[-1][0],max(m[-1][1],b))
            else: m.append((a,b))
        merged[gid]=m
    lanes=sorted(merged.keys(), key=lambda g: min(a for a,_ in merged[g]))
    return merged, lanes, gene_strand, tx_name

def write_tables(tag, chr_, start, end, segs, merged, lane_order, strand_of, name_of):
    peakint={}
    for m in MARK_ORDER:
        fp=f"{BWD}/{m}/{m}_peaks.narrowPeak"
        if not os.path.exists(fp): fp=f"{BWD}/{m}/{m}_peaks.broadPeak"
        inter=[]
        for l in open(fp):
            p=l.split()
            if p[0]!=chr_: continue
            a,b=max(int(p[1]),start),min(int(p[2]),end)
            if a<b: inter.append((a,b))
        peakint[m]=inter
    rows=[["state","state_bp_kb"]+MARK_ORDER]
    for s,e,st in segs:
        row=[f"E{st}",f"{(e-s)/1000:.0f}"]
        for m in MARK_ORDER:
            ov=sum(max(0,min(ee,e)-max(ss,s)) for ss,ee in peakint[m])
            row.append(f"{(ov/1000):.0f}")
        rows.append(row)
    with open(os.path.join(OUTDIR,f"{tag}_state_x_chip.tsv"),"w") as f:
        f.write("\n".join("\t".join(r) for r in rows)+"\n")
    grows=[]
    for gid in lane_order:
        exs=merged[gid]; ex_=expr.get(gid,{})
        grows.append([gid,strand_of[gid],exs[0][0],exs[-1][1],name_of.get(gid,"NA"),
                      ex_.get("exp","NA"),ex_.get("maxfpkm","NA"),ex_.get("meanC","NA"),ex_.get("change","NA")])
    with open(os.path.join(OUTDIR,f"{tag}_genes.tsv"),"w") as f:
        f.write("\t".join(["gene","strand","start","end","protein_name","expressed","max_fpkm","mean_fpkm_Cro","change"])+"\n")
        for r in grows: f.write("\t".join(map(str,r))+"\n")
    print("tables saved:",tag)

TC_MARK_ORDER = ["H3K27ac","H3K4me3","H3K4me1","H3K36me3","H3K27me3","H3K9me3"]

def draw_region_trackc(tag, chr_, start, end, expand, fill):
    """APP-style figure following A485_gene.hic.epi.RNAseq.ipynb: trackc bw/gene/scale tracks
    on a manually stacked layout (trackc 0.0.19 tenon.add cannot stack >1 bottom axes)."""
    if expand: start-=expand*1000; end+=expand*1000
    W=end-start
    segs=load_segments(chr_,start,end)
    merged,lane_order,strand_of,name_of=load_genes(chr_,start,end)

    import pandas as pd
    import trackc as tc
    regions=[f"{chr_}:{start}-{end}"]

    fig=plt.figure(figsize=(11,10))
    AX_X,AX_W=0.13,0.70
    ST_Y,ST_H=0.930,0.032
    CH_H=0.056; STEP=0.068
    mark_y=ST_Y-ST_H-0.012
    ax_state=fig.add_axes([AX_X,ST_Y,AX_W,ST_H])
    ax_bw={}
    for m in TC_MARK_ORDER:
        ax_bw[m]=fig.add_axes([AX_X,mark_y,AX_W,CH_H]); mark_y-=STEP
    ax_rna=fig.add_axes([AX_X,mark_y+STEP-0.004,AX_W,0.058])
    ax_gene=fig.add_axes([AX_X,0.120,AX_W,0.160])
    ax_scale=fig.add_axes([AX_X,0.020,AX_W,0.020])

    for s,e,st in segs:
        ax_state.add_patch(Rectangle((s,0),e-s,1,fc=state_colors[st],ec="none",zorder=3))
    ax_state.set_xlim(start,end); ax_state.set_ylim(0,1)
    ax_state.set_yticks([]); ax_state.set_xticks([])
    ax_state.set_ylabel("ChromHMM",fontsize=10,rotation=0,horizontalalignment="right",
                        verticalalignment="center",labelpad=10)
    ax_state.text(0.004,1.30,f"{chr_}:{start:,}-{end:,}  ({W/1000:.0f} kb)",
                  transform=ax_state.transAxes,fontsize=9)
    leg="\n".join(f"{state_colors[k]} E{k} {state_names[k]}" for k in range(1,11))
    pos=ax_state.get_position()
    fig.text(pos.x1+0.004,(pos.y0+pos.y1)/2,leg,fontsize=8,va="center",ha="left",
             family="monospace")

    bws={m:f"{BWD}/{m}/{m}_FE.bw" for m in TC_MARK_ORDER}
    for m in TC_MARK_ORDER:
        tc.pl.bw_track(bws[m],ax_bw[m],regions=regions,binsize=100,label=m,
                       primary_col="tab:blue")
    tc.pl.bw_track(RNASEQ,ax_rna,regions=regions,binsize=100,label="RNA-seq",
                   vmax=1000.0,primary_col="tab:blue",style="line")

    gb=pd.read_table(BED12,header=None)
    sub=gb[(gb[0]==chr_)&(gb[1]<end)&(gb[2]>start)].copy()
    labels=[str(n) for n in sub[3].tolist()]
    tc.pl.gene_track(ax=ax_gene,bed12=sub,regions=regions,line=1,label="Genes",
                     show_label=labels,ax_on=True,gene_fontszie=9)
    tc.pl.multi_scale_track(ax_scale,regions=regions,scale_adjust="Mb",
                            intervals=1,tick_fontsize=10)

    pdf=os.path.join(OUTDIR,f"{tag}_state_trackc.pdf")
    png=os.path.join(OUTDIR,f"{tag}_state_trackc.png")
    fig.savefig(pdf); fig.savefig(png,dpi=220)
    plt.close(fig)
    print("saved:",pdf); print("saved:",png)
    write_tables(tag,chr_,start,end,segs,merged,lane_order,strand_of,name_of)

def draw_region(tag, chr_, start, end, expand, fill, use_trackc):
    if use_trackc:
        draw_region_trackc(tag, chr_, start, end, expand, fill); return
    if expand: start-=expand*1000; end+=expand*1000
    W=end-start; NBIN=2000
    segs=load_segments(chr_,start,end)
    merged,lane_order,strand_of,name_of=load_genes(chr_,start,end)

    datx=np.linspace(start,end,NBIN)
    bwdata={}
    for m in MARK_ORDER:
        b=pyBigWig.open(f"{BWD}/{m}/{m}_FE.bw")
        vals=b.stats(chr_,start,end,nBins=NBIN,type="mean"); b.close()
        arr=np.array([v if v is not None else 0.0 for v in vals],dtype=float)
        bwdata[m]=(arr,arr.max())
    b=pyBigWig.open(RNASEQ)
    rna_nbin=W//100
    rv=b.stats(chr_,start,end,nBins=rna_nbin,type="mean"); b.close()
    rna=np.array([v if v is not None else 0.0 for v in rv],dtype=float)
    rna_cap=1000.0
    datx_rna=np.linspace(start,end,len(rna))

    # layout: APP style -> taller state bar, lower/thinner chip+rna
    if expand:
        ST_Y,ST_H=0.885,0.0425
        CH_H,STEP,RNA_H,GENE_Y=0.075,0.088,0.080,0.250
        fig_h=9.6
    else:
        ST_Y,ST_H=0.888,0.045
        CH_H,STEP,RNA_H,GENE_Y=0.100,0.112,0.100,0.222
        fig_h=9.2
    fig=plt.figure(figsize=(9,fig_h))
    ax_state=fig.add_axes([0.10,ST_Y,0.72,ST_H])
    axs_bw={}; y=ST_Y-ST_H-0.030
    for m in MARK_ORDER:
        axs_bw[m]=fig.add_axes([0.10,y,0.72,CH_H]); y-=STEP
    y+=STEP
    ax_rna=fig.add_axes([0.10,y-0.005,0.72,RNA_H])
    ax_gene=fig.add_axes([0.10,GENE_Y,0.72,0.150])
    ax_scale=fig.add_axes([0.10,0.045,0.72,0.015])

    for s,e,st in segs:
        ax_state.add_patch(Rectangle((s,0),e-s,1,fc=state_colors[st],ec="none",zorder=3))
    ax_state.set_xlim(start,end); ax_state.set_ylim(0,1)
    ax_state.set_yticks([]); ax_state.set_xticks([])
    ax_state.set_title(f"{chr_}:{start:,}-{end:,}  ({W/1000:.0f} kb)",fontsize=10,loc="left")
    leg="\n".join(f"{state_colors[k]} E{k} {state_names[k]}" for k in range(1,11))
    ax_state.text(1.005,0.5,leg,transform=ax_state.transAxes,fontsize=7.5,va="center",
                  family="monospace")
    for sp in ax_state.spines.values(): sp.set_visible(False)

    fill_color=fill or None
    for m in MARK_ORDER:
        ax=axs_bw[m]; arr,mx=bwdata[m]
        vmax=max(mx*1.05,0.001)
        c=fill_color if fill_color else marks_colors[m]
        ax.fill_between(datx,arr,0,color=c,alpha=0.85,linewidth=0)
        ax.set_xlim(start,end); ax.set_ylim(0,1.08*vmax)
        ax.set_yticks([0,vmax]); ax.set_yticklabels(["0",f"{int(round(vmax))}"],fontsize=6)
        ax.set_xticks([])
        ax.spines["top"].set_visible(False); ax.spines["bottom"].set_visible(False)
        for s in ["left","right"]: ax.spines[s].set_linewidth(0.6)
        ax.set_ylabel(m,fontsize=8,rotation=0,labelpad=30,va="center")
        ax.text(1.005,0.5,f"max {mx:.1f}",transform=ax.transAxes,fontsize=6,va="center")

    vmax=rna_cap
    c=fill_color if fill_color else "#377EB8"
    ax_rna.fill_between(datx_rna,rna,0,color=c,alpha=0.85,linewidth=0)
    ax_rna.set_xlim(start,end); ax_rna.set_ylim(0,vmax)
    ax_rna.set_yticks([0,vmax]); ax_rna.set_yticklabels(["0",f"{int(round(vmax))}"],fontsize=6)
    ax_rna.set_xticks([])
    ax_rna.set_ylabel("Cro RNA-seq\n(RPKM)",fontsize=8,rotation=0,labelpad=30,va="center")
    ax_rna.text(1.005,0.5,f"cap {vmax:.0f}",transform=ax_rna.transAxes,fontsize=6,va="center")
    ax_rna.spines["top"].set_visible(False); ax_rna.spines["bottom"].set_visible(False)

    # genes + scale (custom style for ANXA2/EZR)
    nlane=max(len(lane_order),1)
    for i,gid in enumerate(lane_order):
        exs=merged[gid]; strand=strand_of[gid]
        s=exs[0][0]; e=exs[-1][1]
        ly=1.0-0.15-i*(0.70/nlane)
        color="#1f77b4" if strand=="+" else "#d62728"
        ax_gene.plot([s,e],[ly,ly],color=color,lw=1.5)
        hx=e if strand=="+" else s
        ax_gene.annotate("",xy=(hx,ly+0.10),xytext=(hx,ly-0.10),
                         arrowprops=dict(arrowstyle="-|>",color=color,lw=1.5))
        for xs,xe in exs:
            ax_gene.add_patch(Rectangle((xs,ly-0.07),xe-xs,0.14,fc=color,ec="none"))
        nm=name_of.get(gid,"").replace("_"," ")
        if nm: nm=nm[:28]+"…" if len(nm)>28 else nm
        lab=f"{gid}"+(f" ({nm})" if nm else "")
        ax_gene.text(s,ly-0.32,lab,fontsize=6.5,va="top")
    ax_gene.set_xlim(start,end); ax_gene.set_ylim(0,1)
    ax_gene.set_yticks([]); ax_gene.set_xticks([])
    ax_gene.set_ylabel("Genes",fontsize=8,rotation=0,labelpad=28,va="center")
    for sp in ax_gene.spines.values(): sp.set_visible(False)

    for tick in range(int(math.ceil(start/100000)*100000), end+1, 100000):
        ax_scale.plot([tick,tick],[0,0.4],color="black",lw=0.7)
        if tick>=start: ax_scale.text(tick,0.7,f"{(tick-start)/1e6:.1f}",fontsize=7,ha="center")
    ax_scale.set_xlim(start,end); ax_scale.set_ylim(-0.2,1.4)
    ax_scale.set_yticks([]); ax_scale.set_xticks([])
    ax_scale.set_xlabel(f"{chr_} position (Mb, relative to {start:,})",fontsize=8)
    for sp in ax_scale.spines.values(): sp.set_visible(False)

    pdf=os.path.join(OUTDIR,f"{tag}_state_trackc.pdf")
    png=os.path.join(OUTDIR,f"{tag}_state_trackc.png")
    fig.savefig(pdf); fig.savefig(png,dpi=220)
    plt.close(fig)
    print("saved:",pdf); print("saved:",png)
    write_tables(tag,chr_,start,end,segs,merged,lane_order,strand_of,name_of)

if __name__=="__main__":
    for tag,chr_,start,end,expand,fill,use_trackc in REGIONS:
        draw_region(tag,chr_,start,end,expand,fill,use_trackc)