#!/REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin/python
"""Topic B, Part 1: 12 Nov 2026 post-hoc CI analysis, NOT pre-specified.

Run with /REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin/python -B scripts/posthoc_ci.py.
Only this script and the four explicitly authorised output files are created.
Original analysis scripts are read, never executed wholesale or modified. Python
estimands below replicate the indicated original definitions; the E2 R code
evaluates only the original read/analysis expressions, excluding every write.
No BayesPrism fitting is performed. All numerical libraries use one thread;
the E2 bootstrap uses at most eight single-threaded fork workers.
"""
import os

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
    os.environ[_key] = "1"
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import ast
import hashlib
import json
import subprocess
import sys
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
os.chdir(ROOT)
SEED, NBOOT = 20261112, 10000
DATE = "12 Nov 2026"
ZONE = ZoneInfo("Asia/Shanghai")
LABEL = "Post-hoc 95% percentile bootstrap CIs for the 12 Nov 2026 analysis; not pre-specified"
OUT = ROOT / "results/posthoc_ci"
TARGETS = {
    ROOT / "scripts/posthoc_ci.py", OUT / "posthoc_ci.json", OUT / "posthoc_ci.csv",
    OUT / "POSTHOC_CI_SUMMARY.md", ROOT / "manuscript/drafts/skeleton_v1b_resolution.md",
}
INPUTS = set()
ROWS = []
DECISIONS = [
    ("C01", "Main Fig 1A", "Accounting audit; no CI applies."),
    ("C02", "Main-text A1 detection claim (Fig 2C context)", "Compute paired normal-minus-tumour detection proportion difference."),
    ("C03", "Main Fig 3C and community headline", "Compute no-HP and with-HP rCLR PERMANOVA R2 in the three reported SILVA cohorts."),
    ("C04", "Supporting E1c context; absent from main panels", "CI not computed; retain original effects, P and q only. No main-figure/headline association claim."),
    ("C05", "Supplementary Fig S1A", "Keep existing Methods timing caveat; no change authorised in Part 1."),
    ("C06", "Main Fig 5E; frozen sets in main-text E2a claim", "Compute NES for the 20 plotted Hallmark sets and four frozen sets; discovery counts have no common effect CI."),
    ("C07", "Supporting E2b context; no plotted correlation", "CI not computed; retain effects, P and q only. Programme-score correlations support the supplementary E2 material."),
    ("C08", "Main Fig 5D", "Compute all eight plotted paired M3-gene correlations. Tumour-only MUC5AC is supporting context: CI not computed."),
    ("C09", "Supporting E3c context; supplementary Fig S4A", "CI not computed; retain effects, P and q only. Distinct from the plotted paired H-M estimate."),
    ("C10", "Main-text H-M robustness claim (Fig 5A-C context)", "Compute all three Korean detected-normal subset correlations, as explicitly approved."),
    ("C11", "Main-text audit direction claim (Fig 6A context)", "Compute rank-biserial effects for all 13 subsets and the ten reported single-patient exclusions. Figure intervals themselves concern mixed-model coefficients."),
    ("C12", "Supporting audit diagnostics, not plotted", "CI not computed; retain classifier effect and available breadth P. Classifier P remains unavailable; do not invent it or a breadth effect contrast."),
    ("C13", "Main-text measurement correspondence (Fig 6D context)", "Compute the original overall A1 pooled-sample Spearman correlation using patient-cluster resampling, including unpaired tumours."),
    ("C14", "Supplementary Fig S6B; main Fig 7 interpretability", "Recover deterministic BayesPrism exclusions and check against stored mixture gene names; no refitting."),
    ("C15", "Main Fig 7B", "Compute all nine estimable E9c paired median differences; retain four non-estimable genes without CIs."),
    ("C16", "Main Fig 7C", "Compute all nine estimable E9d correlations; preserve the original natural-log HP change metric."),
    ("C17", "Main Fig 8A A4 association", "Compute the primary and all three registered secondary correlations, using the corrected ssGSEA secondary."),
    ("C18", "Main-text histological-trend headline; plot in Fig S8A", "Compute primary and reported sensitivity JT/Spearman CIs, as explicitly approved despite supplementary plotting."),
]


def path(name):
    p = Path(name)
    if not p.is_absolute():
        p = ROOT / p
    INPUTS.add(p)
    return p


def read_json(name):
    return json.loads(path(name).read_text())


def read_csv(name, **kwargs):
    return pd.read_csv(path(name), **kwargs)


def progress(msg):
    print(msg, flush=True)


def original_functions(name, names):
    """Import only requested pure function ASTs, with no top-level file writes."""
    p = path(name)
    tree = ast.parse(p.read_text())
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in nodes} == set(names)
    namespace = dict(np=np, pd=pd, spearmanr=spearmanr)
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(p), "exec"), namespace)
    return namespace


def add(check, analysis, estimate, original, boots, n, estimand, source, code,
        inputs, *, p=None, q=None, note="", unit="patient", null=0.0,
        affirmative=False, tol=1e-10, extra=None):
    """Reject non-reproduced estimates before accepting any interval."""
    estimate, original = float(estimate), float(original)
    error = abs(estimate - original)
    if not np.isfinite(estimate) or error > tol:
        raise AssertionError(f"Point reproduction failed: {analysis}: {estimate} vs {original}")
    values = np.asarray(boots, dtype=float)
    assert len(values) == NBOOT
    valid = values[np.isfinite(values)]
    if len(valid) < .99 * NBOOT:
        raise AssertionError(f"Too many invalid bootstrap estimates: {analysis}")
    lo, hi = np.percentile(valid, [2.5, 97.5])
    for f in [source.split("::")[0], code, *inputs]:
        path(f)
    row = dict(check=check, analysis=analysis, estimand=estimand, n=int(n),
               estimate=estimate, original_estimate=original, absolute_reproduction_error=error,
               reproduction="PASS", reproduction_tolerance=tol, ci_low=float(lo), ci_high=float(hi),
               confidence_level=.95, method="percentile bootstrap", seed=SEED,
               n_resamples=NBOOT, n_valid=len(valid), n_invalid=NBOOT-len(valid),
               resampling_unit=unit, null_value=null,
               includes_null=bool(lo <= null <= hi), affirmative_claim=affirmative,
               original_P=None if p is None else float(p), original_q=None if q is None else float(q),
               source=source, original_code=code, inputs=inputs, note=note,
               status_label=LABEL, approved_analysis_date="2026-11-12",
               computed_date=datetime.now(ZONE).date().isoformat())
    if extra:
        row.update(extra)
    ROWS.append(row)
    progress(f"{check} {analysis}: {estimate:.8g} [{lo:.8g}, {hi:.8g}]; reproduced")


def indices(n):
    return np.random.default_rng(SEED).integers(0, n, (NBOOT, n))


def rho_rows(x, y):
    x = rankdata(x, axis=-1)
    y = rankdata(y, axis=-1)
    x -= x.mean(axis=-1, keepdims=True)
    y -= y.mean(axis=-1, keepdims=True)
    return (x*y).sum(axis=-1) / np.sqrt((x*x).sum(axis=-1)*(y*y).sum(axis=-1))


def boot_rho(x, y, ix=None):
    x, y = np.asarray(x, float), np.asarray(y, float)
    assert np.isfinite(x).all() and np.isfinite(y).all()
    ix = indices(len(x)) if ix is None else ix
    return np.concatenate([rho_rows(x[v], y[v]) for v in np.array_split(ix, 20)])


def add_rho(check, name, x, y, original, source, code, inputs, **kwargs):
    x, y = np.asarray(x, float), np.asarray(y, float)
    point = spearmanr(x, y).statistic
    # Reproduction is checked before resampling, as well as at serialization.
    assert abs(point - original) < 1e-10, (name, point, original)
    add(check, name, point, original, boot_rho(x, y), len(x), "Spearman rho",
        source, code, inputs, **kwargs)


def accounting():
    manifest = read_csv("data/external/korea16s/ena_full.tsv", sep="\t")
    m = manifest[(manifest.instrument_model == "Illumina MiSeq") & (manifest.library_strategy == "WGS")].copy()
    m["patient"] = m.library_name.str[:-2]
    m["tissue"] = m.library_name.str[-1]
    k = read_csv("results/core/korea_rna_per_sample.csv")
    s = read_csv("results/korea16s/per_sample_silva.csv", index_col=0)
    runs = read_csv("data/external/korea16s/runs16s.tsv", sep="\t", header=None)
    assert set(runs[0]) == set(m.run_accession) == set(s.index)
    assert m.groupby("patient").tissue.agg(set).eq({"N", "T"}).all()
    assert set(m.patient) == set(s.pid)
    overlap = sorted(set(m.patient) & set(k.pid))
    only = sorted(set(m.patient) - set(k.pid))
    # Reconstruct an exact way to obtain the erroneous count. This is evidence
    # of compatibility with a naming-pattern bug, not recovered producer code.
    reconstructed = m.library_name.str.replace(r"_[NT]$", "", regex=True)
    saved = read_json("results/sensitivity/dedup_failbias.json")["korea_16s_vs_rna_overlap"]
    assert reconstructed.nunique() == saved["patients16s"] == 206
    path("scripts/k16_analysis.py")
    result = dict(runs16s=len(m), patients16s=int(m.patient.nunique()), complete_pairs=192,
                  normal_only=0, tumour_only=0, overlap_patients=len(overlap),
                  rna_runs=len(k), rna_patients=int(k.pid.nunique()), rna_complete_pairs=174,
                  rna_tumour_only=sorted(set(k[k.type == "T"].pid)-set(k[k.type == "N"].pid)),
                  sixteen_s_only_patients=only, sixteen_s_only_runs=int(m.patient.isin(only).sum()),
                  original_saved=saved,
                  exact_206_reconstruction={"shared_patient_ids":len(overlap),
                    "parser": "library_name.str.replace(r'_[NT]$', '', regex=True)",
                    "underscore_suffix_runs":int(m.library_name.str.contains(r"_[NT]$").sum()),
                    "hyphen_suffix_runs":int(m.library_name.str.contains(r"-[NT]$").sum()),
                    "unmatched_tissue_specific_library_ids":sorted(m.loc[m.patient.isin(only), "library_name"]),
                    "count":int(reconstructed.nunique())},
                  correct_parser="scripts/k16_analysis.py: runs['pid'] = runs.lib.str[:-2]",
                  producer_provenance="Producer of dedup_failbias.json not present in retained project scripts/logs; referenced sens_v06_s4.py writes sensitivity_s4_results.json only. Available session records contain later reads, not the original write. Removing only underscore tissue suffixes exactly reconstructs 206 because GC13-GC26 use hyphens; this is not recovered historical source code.",
                  conclusion="384 16S runs are 192 complete pairs. The saved 206 field is erroneous; no extra single-tissue 16S patients exist. 178 overlap with RNA; the 14 GC13-GC26 patients account for 28 non-overlapping 16S runs.")
    low=s.groupby("pid").reads.min()<5000
    result["E1_depth_excluded_pairs"]=sorted(low.index[low])
    result["E1_pairs_after_ge5000_both_tissues"]=int((~low).sum())
    assert result["E1_pairs_after_ge5000_both_tissues"]==188
    return result


def a1(k):
    old = read_json("results/core/A1_confirmatory.json")
    w = k.pivot_table(index="pid", columns="type", values="hp_log1p_rpm").dropna()
    det = k.pivot_table(index="pid", columns="type", values="hp_k2_confirmed").loc[w.index] >= 10
    n, t = det.N.to_numpy(), det["T"].to_numpy()
    table = [[int((n&t).sum()),int((n&~t).sum())], [int((~n&t).sum()),int((~n&~t).sum())]]
    assert table == old["mcnemar_detect_ge10"]["table_NposTpos_NposTneg_NnegTpos_NnegTneg"]
    diff = n.astype(float)-t.astype(float)
    original = (table[0][1]-table[1][0])/len(w)
    add("C02", "A1_detection_N_minus_T", diff.mean(), original, diff[indices(len(w))].mean(1), len(w),
        "paired detection proportion difference N-T (threshold >=10 fragments)",
        "results/core/A1_confirmatory.json::mcnemar_detect_ge10", "scripts/core_a1_stats.py",
        ["results/core/korea_rna_per_sample.csv"], p=old["mcnemar_detect_ge10"]["P"],
        note="Original table and both detection proportions reproduced; difference is derived from that stored table, not a previously stored scalar. CI in proportions; multiply by 100 for percentage points.",
        affirmative=True, extra={"table_reproduced":table,"normal_detected":int(n.sum()),"tumour_detected":int(t.sum())})
    d = k[["pid", "s16_clr", "hp_log1p_rpm"]].dropna().reset_index(drop=True)
    original = old["spearman_16S_CLR_vs_log1pRPM"]
    point = spearmanr(d.s16_clr, d.hp_log1p_rpm).statistic
    assert abs(point-original["rho"]) < 1e-10 and len(d)==original["n"]
    clusters = [np.asarray(v) for v in d.groupby("pid",sort=True).indices.values()]
    rng = np.random.default_rng(SEED)
    boots = []
    for _ in range(NBOOT):
        take = np.concatenate([clusters[j] for j in rng.integers(0,len(clusters),len(clusters))])
        boots.append(spearmanr(d.s16_clr.to_numpy()[take],d.hp_log1p_rpm.to_numpy()[take]).statistic)
    add("C13", "A1_overall_RNA_16S", point, original["rho"], boots, len(d), "pooled-sample Spearman rho",
        "results/core/A1_confirmatory.json::spearman_16S_CLR_vs_log1pRPM", "scripts/core_a1_stats.py",
        ["results/core/korea_rna_per_sample.csv"],p=original["P"],unit="patient cluster, all observed tissue rows kept together",
        affirmative=True,extra={"n_patient_clusters":len(clusters),"cluster_size_counts":pd.Series([len(x) for x in clusters]).value_counts().to_dict()})


def community():
    old=read_json("results/step10/E1/E1_results.json")
    for cohort in ["f2","korea16s","yu310127"]:
        gf=f"results/{cohort}/genus_silva.tsv"; mf=f"results/{cohort}/per_sample_silva.csv"
        g=read_csv(gf,sep="\t",index_col=0); m=read_csv(mf,index_col=0)
        if "type" not in m: m["type"]=m.type_bs
        m.pid=m.pid.astype(str); m=m[m.index.isin(g.index)].copy()
        m["total"]=g.loc[m.index].sum(axis=1)
        keep=[]
        for _,grp in m.groupby("pid",sort=True):
            if len(grp)==2 and set(grp.type)=={"N","T"} and grp.total.min()>=5000: keep.extend(grp.index)
        m=m.loc[keep].sort_values(["pid","type"]);g=g.loc[m.index]
        n=len(m)//2; assert n==old[cohort+"_silva"]["n_pairs"]
        # Equivalent sums-of-squares to vegan::adonis2(dist(rclr(g)) ~ type).
        # Paired bootstrap duplicates both tissues; the original prevalence
        # filter and rowwise rCLR centering are recomputed on each resample.
        for drop,label in [(True,"noHP_rCLR"),(False,"withHP_rCLR")]:
            a=g.loc[:,~g.columns.str.contains("g__Helicobacter")] if drop else g
            a=a.to_numpy(float).reshape(n,2,-1)
            with np.errstate(divide="ignore"):
                logs=np.log(a)
            present=a>0; logs=np.where(present,logs,0)
            def statistic(weights):
                mask=np.einsum("i,ijk->k",weights,present)>=.1*(2*n)
                l=logs[:,:,mask]; hit=present[:,:,mask]
                mu=l.sum(axis=2)/hit.sum(axis=2)
                x=np.where(hit,l-mu[:,:,None],0)
                means=np.einsum("i,ijk->jk",weights,x)/n
                overall=means.mean(axis=0)
                ss=np.einsum("i,ijk->",weights,(x-overall)**2)
                return float((n*((means-overall)**2).sum())/ss)
            original=old[cohort+"_silva"]["E1a"][label]
            point=statistic(np.ones(n));assert abs(point-original["R2"])<1e-10
            rng=np.random.default_rng(SEED)
            boots=[statistic(np.bincount(rng.integers(0,n,n),minlength=n)) for _ in range(NBOOT)]
            add("C03",f"E1a_{cohort}_{label}",point,original["R2"],boots,n,"PERMANOVA R2 (Euclidean rCLR)",
                f"results/step10/E1/E1_results.json::{cohort}_silva.E1a.{label}","scripts/s10_e1.R",[gf,mf],p=original["P"],
                note="Original pair/depth rule, HP handling, >=10% prevalence and rCLR definition replicated. R2 is nonnegative; its percentile CI excluding zero does not replace the permutation test or resolve the betadisper caveat.")


def e3_hm(k):
    old=read_json("results/step10/E3/E3_results.json")
    e=read_csv("results/core/korea_gene_log2tpm.csv.gz",index_col=0).T
    hm=read_csv("results/core/hm_deltas_korea.csv",dtype={"pid":str}).set_index("pid")
    w=k.pivot_table(index="pid",columns="type",values="hp_log10rpm1").dropna()
    ids=list(w.index);nn=[p+"_N" for p in ids];tt=[p+"_T" for p in ids]
    hp=w["T"]-w.N
    assert np.allclose(hp,hm.loc[ids,"dHP"])
    # Original E3d uses log2 bulk expression differences and H-M log10 HP.
    for gene in ["MUC5AC","TFF1","TFF2","MUC6","MUC1","FUT2","FUT3","B3GNT6"]:
        d=e.loc[tt,gene].to_numpy()-e.loc[nn,gene].to_numpy();o=old["E3d"][gene]
        q=next(x["q"] for x in old["BH_family_E3"] if x["item"]=="E3d" and x["test"]==gene)
        add_rho("C08",f"E3d_{gene}",hp,d,o["rho"],f"results/step10/E3/E3_results.json::E3d.{gene}",
            "scripts/s10_e3.py",["results/core/korea_gene_log2tpm.csv.gz","results/core/korea_rna_per_sample.csv"],p=o["P"],q=q,
            affirmative=gene in ["TFF1","TFF2"])
    sens=read_json("results/sensitivity/sensitivity_s4_results.json")["hm_v05_sensitivities"]["korea"]
    # Rebuild the original frozen marker scores, then match the saved deltas.
    f=original_functions("scripts/sens_v06_s4.py",["zm","scores"])
    tree=ast.parse(path("scripts/sens_v06_s4.py").read_text())
    for node in tree.body:
        if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ["IM","GD","M2","M3","EPI","STR","IMM"] for t in node.targets):
            exec(compile(ast.Module(body=[node],type_ignores=[]),"scores constants","exec"),f)
    scores=f["scores"](e)
    normal_counts=k.set_index("lib").loc[nn,"hp_k2_confirmed"].to_numpy()
    keep=normal_counts>=10
    for name in ["M1","M2","M3"]:
        d=scores.loc[tt,name].to_numpy()-scores.loc[nn,name].to_numpy()
        assert np.allclose(d,hm.loc[ids,"d"+name])
        o=sens["sens3_"+name]
        add_rho("C10",f"HM_detected_normal_{name}",hp.to_numpy()[keep],d[keep],o["rho"],
            f"results/sensitivity/sensitivity_s4_results.json::hm_v05_sensitivities.korea.sens3_{name}",
            "scripts/sens_v06_s4.py",["results/core/korea_gene_log2tpm.csv.gz","results/core/korea_rna_per_sample.csv","results/core/hm_deltas_korea.csv"],
            p=o["P"],affirmative=name in ["M1","M2"],note="Scores rebuilt and matched to saved H-M deltas. Patient bootstrap conditional on the original cohort-standardized scores and >=10 normal-fragment subset; scores are not reselected.")


def audit():
    d=read_csv("results/audit/audit_per_run.csv",dtype={"pid":str})
    old=read_json("results/audit/audit_results.json")
    w=d.pivot_table(index="pid",columns="type",values="hp_log1p_rpm").dropna()
    top=(w.N-w["T"]).abs().sort_values(ascending=False).index
    mx=d.groupby("pid").hp_k2_confirmed.max().sort_values(ascending=False).index
    sources=[]
    for j,o in enumerate(old["subsets"]):
        dx=d.copy(); col=o["count_col"]
        if col=="hp_k2_restr": dx[col]=dx.hp_k2_confirmed.where(dx.k2out_available==True)
        if j==10: dx=dx[~dx.pid.isin(top[:5])]
        if j==11: dx=dx[~dx.pid.isin(top[:10])]
        if j==12: dx=dx[~dx.pid.isin(mx[:10])]
        sources.append((f"audit_subset_{j:02d}",dx,col,o["rank_biserial"],o["P_one_sided_N_gt_T"],o["subset"],f"subsets[{j}]"))
    for o in old["loo_single_top10"]["rows"]:
        sources.append(("audit_drop_"+o["pid"],d[d.pid!=o["pid"]],"hp_k2_confirmed",o["rb"],o["P1"],"single-patient exclusion", "loo_single_top10.rows:"+o["pid"]))
    for name,dx,col,orig,p,note,key in sources:
        dx=dx.copy();dx["y"]=np.log1p(dx[col]/dx.salmon_mapped*1e6)
        v=dx.pivot_table(index="pid",columns="type",values="y").dropna();delta=(v.N-v["T"]).to_numpy()
        def rb(x):
            ranks=rankdata(np.where(x==0,np.nan,abs(x)),axis=-1,nan_policy="omit")
            return np.nansum(ranks*np.sign(x),axis=-1)/np.nansum(ranks,axis=-1)
        point=rb(delta);assert abs(point-orig)<1e-10
        boots=np.concatenate([rb(delta[ix]) for ix in np.array_split(indices(len(delta)),20)])
        add("C11",name,point,orig,boots,len(delta),"paired rank-biserial r, N-T",
            "results/audit/audit_results.json::"+key,"scripts/audit_stage2.py",["results/audit/audit_per_run.csv"],p=p,
            note=note+"; resample all eligible patient pairs including zero differences, then omit zeros only when calculating signed ranks. Original exclusion set held fixed.",affirmative=True)


def e9(k):
    old=read_json("results/step10/E9/E9_results.json")["genes"]
    f=original_functions("scripts/s10_e7_e9_analysis.py",["epi_cpm"])
    ep=f["epi_cpm"](path("results/step10/E7/korea_epithelial_exp.csv.gz"))
    w=k.pivot_table(index="pid",columns="type",values="hp_log1p_rpm").dropna()
    ids=list(w.index);tt=[p+"_T" for p in ids];nn=[p+"_N" for p in ids]
    for gene,o in old.items():
        if gene.startswith("_") or not isinstance(o["E9c"],dict):continue
        delta=ep.loc[tt,gene].to_numpy()-ep.loc[nn,gene].to_numpy()
        point=np.median(delta);assert abs(point-o["E9c"]["median_diff"])<1e-10
        add("C15",f"E9c_{gene}",point,o["E9c"]["median_diff"],np.median(delta[indices(len(delta))],axis=1),len(delta),
            "median paired T-N log2(epithelial CPM+1)",f"results/step10/E9/E9_results.json::genes.{gene}.E9c",
            "scripts/s10_e7_e9_analysis.py",["results/step10/E7/korea_epithelial_exp.csv.gz","results/core/korea_rna_per_sample.csv"],
            p=o["E9c"]["P"],q=o["E9c"]["q_BH"],affirmative=gene not in ["MUC1","FUT2"],
            note="Conditional on stored BayesPrism epithelial expression; no BayesPrism refit. Same epithelial CPM denominator as original epi_cpm().")
        add_rho("C16",f"E9d_{gene}",w["T"]-w.N,delta,o["E9d"]["rho"],f"results/step10/E9/E9_results.json::genes.{gene}.E9d",
            "scripts/s10_e7_e9_analysis.py",["results/step10/E7/korea_epithelial_exp.csv.gz","results/core/korea_rna_per_sample.csv"],
            p=o["E9d"]["P"],q=o["E9d"]["q_BH"],affirmative=gene in ["TFF1","TFF2","FUT2","A4GNT","CEACAM6"],
            note="HP predictor is T-N natural log1p(RPM), not the H-M log10 metric. Conditional on stored epithelial expression.")


def a4():
    d=read_csv("results/step10/A4full/korea_normal_scores_full.csv",index_col=0)
    old=read_json("results/step10/A4full/A4_full_results.json")
    s=read_csv("results/step10/A4full/IM_ssgsea_normals.csv").set_index("lib")
    assert np.allclose(d.IM_ssgsea,s.loc[d.index,"IM_ssgsea"])
    for o,name in [(old["primary"],"primary_CLR_IM")]+[(x,"secondary_"+x["x"]+"_"+x["y"]) for x in old["secondary"]]:
        v=d[[o["x"],o["y"]]].dropna()
        key="primary" if name.startswith("primary") else "secondary[y="+o["y"]+",x="+o["x"]+"]"
        add_rho("C17","A4_"+name,v.iloc[:,0],v.iloc[:,1],o["rho"],
            "results/step10/A4full/A4_full_results.json::"+key,"scripts/s10_a4full_e10.py",
            ["results/step10/A4full/korea_normal_scores_full.csv","results/step10/A4full/IM_ssgsea_normals.csv","scripts/s10_a4full_ssgsea.R"],
            p=o["P"],q=o.get("q_BH"),affirmative=o["x"]!="dCLR",
            note="Original stored, cohort-standardized normal scores; corrected GSVA ssGSEA secondary, not the superseded mean-rank score. Conditional on original score construction.")


def e12():
    d=read_csv("results/prjeb9749/E12_E13/per_sample.csv");d=d[d.QC_pass].copy()
    old=read_json("results/prjeb9749/E12_E13/E12_E13_results.json")
    jt=original_functions("scripts/e12_e13_prjeb9749.py",["jt"])["jt"]
    variants=[("primary",d,"y","group5_order",old["E12_primary"]),
        ("S1_factor4",d,"y","factor4_order",old["E12_sensitivity"]["S1_official_factor4_JT"]),
        ("S3_fastp_denominator",d,"y_k2only","group5_order",old["E12_sensitivity"]["S3_kraken2_only_JT"]),
        ("S4_drop_EA1",d[d.run!="ERR950183"],"y","group5_order",old["E12_sensitivity"]["S4_drop_EA1_ERR950183_JT"])]
    for name,v,col,gcol,o in variants:
        key={"primary":"E12_primary","S1_factor4":"E12_sensitivity.S1_official_factor4_JT",
             "S3_fastp_denominator":"E12_sensitivity.S3_kraken2_only_JT","S4_drop_EA1":"E12_sensitivity.S4_drop_EA1_ERR950183_JT"}[name]
        y=v[col].to_numpy();g=v[gcol].to_numpy();groups=[np.flatnonzero(g==level) for level in np.unique(g)]
        j=jt(y,g);assert j==o["J_obs"]
        jmax=sum(len(a)*len(b) for z,a in enumerate(groups) for b in groups[z+1:]);assert jmax==o["J_max"]
        ix=indices(len(v))
        boots=np.array([jt(y[i],g[i]) for i in ix])
        boot_max=np.array([sum(np.sum(g[i]==a)*np.sum(g[i]==b) for z,a in enumerate(np.unique(g)) for b in np.unique(g)[z+1:]) for i in ix])
        add("C18","E12_JT_"+name,j/jmax,o["J_over_Jmax"],boots/boot_max,len(v),"JT J/Jmax, decreasing direction",
            "results/prjeb9749/E12_E13/E12_E13_results.json::"+key,"scripts/e12_e13_prjeb9749.py",["results/prjeb9749/E12_E13/per_sample.csv"],
            p=o["P_one_sided_decreasing"],null=.5,unit="subject (ordinary bootstrap, observed group and outcome resampled together)",
            note="Original JT function imported without top-level execution. Each bootstrap recalculates Jmax for its resampled group sizes. Null is 0.5 (not zero); no new permutation P computed. Raw J CI is also provided but depends on resampled group sizes.",
            extra={"J_observed_reproduced":j,"J_max_reproduced":jmax,"J_ci_low":float(np.percentile(boots,2.5)),"J_ci_high":float(np.percentile(boots,97.5))})
        if name in ["primary","S1_factor4"]:
            sp=o["spearman_effect"] if name=="primary" else old["E12_sensitivity"]["S2_spearman_ordinal"]["factor4"]
            point=spearmanr(y,g).statistic;assert abs(point-sp["rho"])<1e-10
            add("C18","E12_Spearman_"+name,point,sp["rho"],boot_rho(y,g,ix),len(v),"Spearman rho, increasing histological order",
                "results/prjeb9749/E12_E13/E12_E13_results.json::"+("E12_primary.spearman_effect" if name=="primary" else "E12_sensitivity.S2_spearman_ordinal.factor4"),"scripts/e12_e13_prjeb9749.py",["results/prjeb9749/E12_E13/per_sample.csv"],
                p=sp["P_two"],unit="subject (ordinary bootstrap, observed group and outcome resampled together)",
                note="Five-group S2 repeats the primary descriptive rho and shares its CI; original one-sided sensitivity P is retained in the source. No separate duplicate resampling.")


def run_r(code, on_chunk=None):
    """R receives code through stdin, preventing auxiliary source/output files."""
    proc=subprocess.Popen(["/REPLACE_WITH_ENVIRONMENT_ROOT/envs/rbio/bin/Rscript","--vanilla","-"],
                          stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
    proc.stdin.write(code);proc.stdin.close()
    result=None
    for line in proc.stdout:
        if line.startswith("__POSTHOC_BATCH__"):
            if on_chunk is not None:on_chunk(json.loads(line[len("__POSTHOC_BATCH__"):]))
        elif line.startswith("__POSTHOC_JSON__"):
            result=json.loads(line[len("__POSTHOC_JSON__"):])
        else:
            progress("R: "+line.rstrip())
    if proc.wait()!=0 or result is None:
        raise RuntimeError("Read-only R calculation failed; see console output")
    return result


BAYES_R = r'''
suppressPackageStartupMessages({library(BayesPrism);library(SummarizedExperiment);library(jsonlite)})
B <- getwd(); marker<-readLines(file.path(B,'results/step10/E7/marker_genes.txt'))
out<-list(package_version=as.character(packageVersion('BayesPrism')),filter_definition=paste(deparse(BayesPrism:::filter.bulk.outlier),collapse='\n'))
check<-function(mix,tag){
 x<-cleanup.genes(input=mix,input.type='count.matrix',species='hs',gene.group=c('Rb','Mrp','other_Rb','chrM','MALAT1','chrX','chrY'),exp.cells=5)
 before_pc<-colnames(x); x<-select.gene.type(x,gene.type='protein_coding')
 y<-BayesPrism:::filter.bulk.outlier(mixture=x,outlier.cut=.01,outlier.fraction=.1)
 exc<-setdiff(colnames(x),colnames(y)); retained<-intersect(marker,colnames(y))
 bp<-readRDS(file.path(B,paste0('results/step10/E7/bp_',tag,'.rds')))
 # Read existing prism mixture, never call new.prism() or run.prism().
 stored<-colnames(bp@prism@mixture)
 stopifnot(setequal(stored,retained),length(exc)==if(tag=='korea')16 else 8)
 lg<-jsonlite::fromJSON(gsub(': NaN',': null',paste(readLines(file.path(B,'results/step10/E7/E7_results.json')),collapse='\n'),fixed=TRUE))$E7b_model2_gene_presence
 score_missing<-if(tag=='korea')lapply(lg,function(z){gs<-z$missing;list(missing=gs,removed_as_outlier=intersect(gs,exc),absent_from_marker_set=setdiff(gs,marker))}) else NULL
 list(n_samples=nrow(mix),n_raw_genes=ncol(mix),n_after_cleanup=length(before_pc),n_after_protein_coding=ncol(x),n_outlier_genes=length(exc),
      outlier_genes=exc,outlier_marker_genes=intersect(exc,marker),n_shared_retained=length(retained),stored_mixture_gene_set_identical=TRUE,
      marker_genes_missing_before_outlier=setdiff(marker,colnames(x)),score_missing_breakdown=score_missing,
      outlier_sample_fractions=setNames(as.list(colMeans(x[,exc,drop=FALSE]/rowSums(x)>.01)),exc))
}
kb<-t(as.matrix(read.csv(file.path(B,'results/step10/E7/input/korea_bulk_counts.csv.gz'),row.names=1,check.names=FALSE)))
out$korea<-check(kb,'korea');rm(kb);invisible(gc())
se<-readRDS(file.path(B,'data/host/star/stad_star_se.rds'));pid<-read.csv(file.path(B,'results/core/hm_deltas_tcga.csv'))$pid
cd<-as.data.frame(colData(se));s16<-substr(cd$barcode,1,16)
sp<-read.csv(file.path(B,'data/external/spectrum2023/stad_helicobacter_per_sample.csv'),row.names=1)
pick<-function(p,ty){ids<-rownames(sp)[substr(rownames(sp),1,12)==p & substr(rownames(sp),14,15)==ty];k<-substr(ids[1],1,16);which(s16==k)[1]}
cols<-c(sapply(pid,pick,ty='11'),sapply(pid,pick,ty='01'));stopifnot(!any(is.na(cols)))
m<-assay(se,'unstranded')[,cols];colnames(m)<-c(paste0(pid,'_N'),paste0(pid,'_T'));m<-rowsum(m,rowData(se)$gene_name)
out$tcga<-check(t(m),'tcga')
cat('__POSTHOC_JSON__',toJSON(out,auto_unbox=TRUE,digits=NA,null='null'), '\n',sep='')
'''


E2_ACCEL = r'''
import numpy as np
from scipy.stats import rankdata

def trim_ranks(a,lo,hi):
    # Exact average-rank interval membership, using selection instead of
    # sorting every gene. Boundary ties are included by their average rank.
    edges=np.partition(a,[int(lo)-1,int(hi)-1])[[int(lo)-1,int(hi)-1]]
    lower,upper=edges
    keep=(a>=lower)&(a<=upper)
    lower_rank=(np.sum(a<lower)+1+np.sum(a<=lower))/2
    upper_rank=(np.sum(a<upper)+1+np.sum(a<=upper))/2
    if lower_rank<lo:keep &= a!=lower
    if upper_rank>hi:keep &= a!=upper
    return keep

def tmm_numpy(counts):
    # Literal edgeR .calcFactorTMM formulas, evaluated once per distinct
    # duplicated sample. Rank trimming retains R's average-rank tie rule.
    x=np.asarray(counts,dtype=float)
    lib=x.sum(axis=0)
    x=x[np.any(x>0,axis=1)]
    f75=np.quantile(x,.75,axis=0)/lib
    refcol=(np.argmax(np.sqrt(x).sum(axis=0)) if np.median(f75)<1e-20
            else np.argmin(abs(f75-f75.mean())))
    ref=x[:,refcol];nr=lib[refcol]
    out=np.ones(x.shape[1]);cache={}
    for j in range(x.shape[1]):
        obs=x[:,j];no=lib[j]
        key=obs.tobytes()
        if key in cache:out[j]=cache[key];continue
        with np.errstate(divide='ignore',invalid='ignore'):
            lr=np.log2((obs/no)/(ref/nr))
            ae=(np.log2(obs/no)+np.log2(ref/nr))/2
            v=(no-obs)/no/obs+(nr-ref)/nr/ref
        good=np.isfinite(lr)&np.isfinite(ae)&(ae>-1e10)
        lr=lr[good];ae=ae[good];v=v[good]
        if np.max(abs(lr))<1e-6:out[j]=1.;cache[key]=1.;continue
        n=len(lr);loL=np.floor(n*.3)+1;hiL=n+1-loL;loS=np.floor(n*.05)+1;hiS=n+1-loS
        keep=trim_ranks(lr,loL,hiL)&trim_ranks(ae,loS,hiS)
        f=np.sum(lr[keep]/v[keep])/np.sum(1/v[keep]);f=0. if np.isnan(f) else f
        out[j]=2**f;cache[key]=out[j]
    return out/np.exp(np.mean(np.log(out)))

def weighted_fit_numpy(expression,weights,design):
    # limma lm.series weighted least squares, with an invertible centering
    # and scaling of design columns for numerical conditioning. Only the
    # HP coefficient / unscaled SE and residual variances enter its t score.
    y=np.asarray(expression);w=np.asarray(weights);X=np.array(design,copy=True)
    X[:,1:]=(X[:,1:]-X[:,1:].mean(axis=0))/X[:,1:].std(axis=0,ddof=1)
    xx=(X[:,:,None]*X[:,None,:]).reshape(len(X),-1)
    cross=(w@xx).reshape(len(y),X.shape[1],X.shape[1])
    rhs=(w*y)@X
    beta=np.linalg.solve(cross,rhs[...,None])[...,0]
    inverse=np.linalg.inv(cross)
    residual=y-beta@X.T
    df=len(X)-X.shape[1]
    variance=(w*residual**2).sum(axis=1)/df
    base=beta[:,1]/np.sqrt(inverse[:,1,1])
    return {'base':base,'variance':variance,'df':float(df)}

def voom_prepare_numpy(counts,norm_factors,design):
    # limma::voom's logCPM, unweighted fit and mean-variance coordinates.
    # R's unchanged lowess and adaptive span are applied to sx/sy below.
    x=np.asarray(counts);X=np.asarray(design);lib=x.sum(axis=0)*np.asarray(norm_factors)
    y=np.log2((x+.5)/(lib+1)*1e6)
    Q,R=np.linalg.qr(X,mode='reduced')
    coef=np.linalg.solve(R,Q.T@y.T).T
    fitted=coef@X.T
    sigma=np.sqrt(((y-fitted)**2).sum(axis=1)/(len(X)-X.shape[1]))
    sx=y.mean(axis=1)+np.log2(lib+1).mean()-np.log2(1e6)
    sy=np.sqrt(sigma)
    logcount=np.log2(1e-6*np.exp2(fitted)*(lib+1))
    return {'expression':y,'sx':sx,'sy':sy,'logcount':logcount}

def voom_finish_numpy(logcount,lx,ly):
    # R approxfun(...,rule=2,ties=list('ordered',mean)) on the lowess curve.
    lx=np.asarray(lx);ly=np.asarray(ly)
    ux,first,count=np.unique(lx,return_index=True,return_counts=True)
    uy=np.add.reduceat(ly,first)/count
    return 1/np.interp(np.asarray(logcount),ux,uy)**4
'''


E2_R = r'''
suppressPackageStartupMessages({library(limma);library(edgeR);library(fgsea);library(jsonlite);library(BiocParallel);library(data.table)})
data.table::setDTthreads(1L)
p<-parse('scripts/s10_e2.R')
# Evaluate original expressions through the t table and gene sets, with every
# write expression omitted; no E2b, TCGA or plotting code is executed.
for(i in c(2:24,26:31)) eval(p[[i]])
stopifnot(length(p)==42,grepl('write.csv',deparse(p[[25]])[1]),grepl('fgsea',deparse(p[[32]])[1]))
nes<-function(st){
 pp<-fgsea:::preparePathwaysAndStats(sets,st,5,500,1,'std')
 st<-pp$stats;gs<-do.call(rbind,lapply(pp$filtered,fgsea::calcGseaStat,stats=st,returnLeadingEdge=TRUE,scoreType='std'))
 seeds<-sample.int(1e9,1)
 z<-fgsea:::fgseaSimpleImpl(pathwayScores=unlist(gs[,'res']),pathwaysSizes=pp$sizes,pathwaysFiltered=pp$filtered,
     leadingEdges=mapply('[',list(names(st)),gs[,'leadingEdge'],SIMPLIFY=FALSE),permPerProc=1000,seeds=seeds,
     toKeepLength=length(pp$filtered),stats=st,BPPARAM=BiocParallel::SerialParam(),scoreType='std')
 # This is the unchanged NES branch inside fgseaMultilevel. The subsequent
 # adaptive P-value branch does not change NES; do not recompute P-values.
 z$NES[ifelse(z$ES>=0,z$nGeZero,z$nLeZero)<10]<-NA_real_
 setNames(z$NES,z$pathway)
}
old<-read.csv('results/step10/E2/E2a_fgsea.csv');oldt<-read.csv('results/step10/E2/E2a_limma_normals.csv.gz',row.names=1)
t_error<-max(abs(tt[rownames(oldt),'t']-oldt$t))
set.seed(20261030);point<-nes(st);nes_error<-max(abs(point[old$pathway]-old$NES))
stopifnot(t_error<1e-8,nes_error<1e-8)
cat('E2 original t/NES reproduced; max errors',t_error,nes_error,'\n');flush.console()
Sys.setenv(RETICULATE_PYTHON='/REPLACE_WITH_ENVIRONMENT_ROOT/envs/shotgun/bin/python',RETICULATE_USE_MANAGED_VENV='no',PYTHONDONTWRITEBYTECODE='1')
# The absolute Python already works without conda activation. Avoid an
# unnecessary mamba cache write; this alters only the current R process.
loadNamespace('reticulate')
assignInNamespace('python_munge_path',function(python){old<-Sys.getenv('PATH');Sys.setenv(PATH=paste(dirname(python),old,sep=.Platform$path.sep));old},ns='reticulate')
reticulate::py_run_string(__PY_CODE__)
positive<-old[grepl('^HALLMARK_',old$pathway)&old$padj<.1&old$NES>0,];positive<-positive[order(-positive$NES),]
negative<-old[grepl('^HALLMARK_',old$pathway)&old$padj<.1&old$NES<0,]
targets<-c(head(positive$pathway,15),negative$pathway,old$pathway[grepl('^V05',old$pathway)])
stopifnot(length(targets)==24)
Kn_original<-Kn
# Each resample repeats original normal-score construction, gene filtering,
# TMM normalization, voom weights, lmFit/eBayes and gene ranking.
set.seed(20261112);draws<-replicate(10000,sample.int(nrow(Kn),nrow(Kn),replace=TRUE),simplify=FALSE)
resample_seeds<-sample.int(1e9,10000)
original_one<-function(j){
 env<-new.env(parent=globalenv());env$Kn<-Kn_original[draws[[j]],,drop=FALSE]
 for(i in 14:24) eval(p[[i]],envir=env)
 setNames(env$tt$t,rownames(env$tt))
}
fast_stats<-function(j){
 env<-new.env(parent=globalenv());env$Kn<-Kn_original[draws[[j]],,drop=FALSE]
 for(i in 14:20) eval(p[[i]],envir=env)
 env$y$samples$norm.factors<-as.numeric(reticulate::py$tmm_numpy(env$y$counts))
 vp<-reticulate::py$voom_prepare_numpy(env$y$counts,env$y$samples$norm.factors,env$des)
 span<-limma::chooseLowessSpan(nrow(env$y),small.n=50,min.span=.3,power=1/3)
 curve<-stats::lowess(as.numeric(vp$sx),as.numeric(vp$sy),f=span)
 weights<-reticulate::py$voom_finish_numpy(vp$logcount,curve$x,curve$y)
 fit<-reticulate::py$weighted_fit_numpy(vp$expression,weights,env$des)
 shrink<-limma::squeezeVar(as.numeric(fit$variance),as.numeric(fit$df))
 setNames(as.numeric(fit$base)/sqrt(shrink$var.post),rownames(env$y$counts))
}
# Validate the accelerated algebra against the original R path on the full
# data and five deterministic bootstrap draws, before producing any CI.
draws[[10001]]<-seq_len(nrow(Kn_original))
validation<-numeric()
for(j in c(10001,1:5)){
 ref<-if(j==10001)st else original_one(j)
 got<-fast_stats(j);err<-max(abs(ref-got[names(ref)]));validation<-c(validation,err)
 stopifnot(identical(names(ref),names(got)),err<1e-8)
 set.seed(20261030);rn<-nes(ref);set.seed(20261030);gn<-nes(got)
 stopifnot(max(abs(rn-gn))<1e-8)
}
cat('E2 accelerated exact-formula validation t/NES PASS; max t error',max(validation),'\n');flush.console()
one<-function(j){
 stj<-fast_stats(j)
 set.seed(resample_seeds[j]);v<-nes(stj)
 v[targets]
}
all<-matrix(NA_real_,10000,length(targets),dimnames=list(NULL,targets))
resume<-jsonlite::fromJSON(__RESUME_JSON__)
done<-0L
if(length(resume$values)){
 stopifnot(identical(resume$targets,targets));done<-nrow(resume$values);all[seq_len(done),]<-resume$values
}
if(done<10000)for(start in seq(done+1,10000,by=80)){
 tick<-proc.time()[3]
 jobs<-start:min(start+79,10000)
 ans<-parallel::mclapply(jobs,one,mc.cores=8,mc.set.seed=FALSE,mc.preschedule=TRUE)
 stopifnot(!any(vapply(ans,inherits,logical(1),'try-error')))
 all[jobs,]<-do.call(rbind,ans)
 cat('__POSTHOC_BATCH__',toJSON(list(start=min(jobs),targets=targets,values=unname(all[jobs,,drop=FALSE])),auto_unbox=TRUE,matrix='rowmajor',digits=NA,na='null'), '\n',sep='');flush.console()
 cat('E2 patient resamples',max(jobs),'/ 10000; last batch seconds',proc.time()[3]-tick,'\n');flush.console()
}
rr<-lapply(targets,function(nm){v<-all[,nm];list(name=nm,estimate=unname(point[nm]),original=old$NES[match(nm,old$pathway)],
 ci=as.numeric(quantile(v,c(.025,.975),na.rm=TRUE,type=7)),n_valid=sum(is.finite(v)),P=old$pval[match(nm,old$pathway)],q=old$padj[match(nm,old$pathway)])})
cat('__POSTHOC_JSON__',toJSON(list(rows=rr,t_max_error=t_error,NES_max_error=nes_error,accelerated_validation_max_t_error=max(validation),n=174,versions=sapply(c('limma','edgeR','fgsea'),function(x)as.character(packageVersion(x)))),auto_unbox=TRUE,digits=NA), '\n',sep='')
'''


def e2(meta):
    files=["scripts/s10_e2.R","results/core/korea_rna_per_sample.csv","results/step10/E7/input/korea_bulk_counts.csv.gz",
        "results/core/korea_gene_log2tpm.csv.gz","results/step10/E2/E2a_fgsea.csv","results/step10/E2/E2a_limma_normals.csv.gz",
        "/REPLACE_WITH_REFERENCE_ROOT/msigdb/h.all.v2025.1.Hs.symbols.gmt"]
    for f in files:path(f)
    state=meta.setdefault("E2_bootstrap_checkpoint",{"targets":[],"values":[]})
    def chunk(o):
        assert o["start"]==len(state["values"])+1
        state["targets"]=o["targets"];state["values"].extend(o["values"])
        write_checkpoint(meta)
    code=E2_R.replace("__PY_CODE__",json.dumps(E2_ACCEL)).replace("__RESUME_JSON__",json.dumps(json.dumps(state)))
    result=run_r(code,on_chunk=chunk)
    for o in result["rows"]:
        assert o["n_valid"]>=9900 and abs(o["estimate"]-o["original"])<1e-8
        lo,hi=o["ci"]
        ROWS.append(dict(check="C06",analysis="E2a_"+o["name"],estimand="fgsea normalized enrichment score (NES)",n=174,
            estimate=o["estimate"],original_estimate=o["original"],absolute_reproduction_error=abs(o["estimate"]-o["original"]),
            reproduction="PASS",reproduction_tolerance=1e-8,ci_low=lo,ci_high=hi,confidence_level=.95,method="percentile bootstrap",
            seed=SEED,n_resamples=NBOOT,n_valid=o["n_valid"],n_invalid=NBOOT-o["n_valid"],resampling_unit="patient (one normal per patient)",
            null_value=0.,includes_null=bool(lo<=0<=hi),affirmative_claim=True,
            affirmative_claim_basis="Significant enrichment displayed in main Fig 5E or named in the main-text E2a claim; original selection was padj <0.10.",
            original_P=o["P"],original_q=o["q"],source="results/step10/E2/E2a_fgsea.csv::"+o["name"],original_code="scripts/s10_e2.R",inputs=files,
            note="Original R expressions plus algebraically identical NumPy TMM, voom matrix algebra and weighted least squares (validated against original R on full data and five bootstrap draws); installed R lowess, limma adaptive span/squeezeVar and fgsea NES code. Original seed used only for reproduction. 10,000 patient resamples, seed 20261112; NES inner null draws have deterministic sub-seeds. Adaptive P-value branch omitted because it does not alter NES. Fixed original plotted/headline sets; no reselection.",
            status_label=LABEL,approved_analysis_date="2026-11-12",computed_date=datetime.now(ZONE).date().isoformat()))
    del meta["E2_bootstrap_checkpoint"]
    return {k:v for k,v in result.items() if k!="rows"}


def snapshot():
    return {str(p.relative_to(ROOT)):(p.stat().st_size,p.stat().st_mtime_ns)
            for p in ROOT.rglob("*") if p.is_file() and p not in TARGETS}


def finite_json(x):
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)


def write_checkpoint(meta):
    OUT.mkdir(exist_ok=True)
    (OUT/"posthoc_ci.json").write_text(json.dumps({**meta,"intervals":ROWS},indent=2,default=finite_json,allow_nan=False)+"\n")


def fmt(x):
    return f"{x:.6g}"


def supporting_results():
    """Copy the effects/P values for analyses deliberately receiving no CI."""
    e1=read_json("results/step10/E1/E1_results.json")
    e2=read_json("results/step10/E2/E2_results.json")
    e3=read_json("results/step10/E3/E3_results.json")
    au=read_json("results/audit/audit_results.json")
    return {
        "C04":dict(status="CI not computed",source="results/step10/E1/E1_results.json",reported_Gemella=e1["f2_silva"]["E1c"]["Gemella"]),
        "C07":dict(status="CI not computed",source="results/step10/E2/E2_results.json",paired=e2["E2b"]["paired_delta"],tumour=e2["E2b"]["tumour"]),
        "C08_tumour_context":dict(status="CI not computed",source="results/step10/E3/E3_results.json",tumour_MUC5AC=e3["E3d"]["Oosterlinck_tumour_crosssectional_HP_vs_MUC5AC"]),
        "C09":dict(status="CI not computed",source="results/step10/E3/E3_results.json",normal_HP_delta_M1=e3["E3c"]["i_normalHP_vs_dM_M1"],normal_HP_IM=e3["E3c"]["ii_withinNormal_HP_vs_IM_A4"]),
        "C12":dict(status="CI not computed; classifier P unavailable",source="results/audit/audit_results.json",
                   classifier_rho=au["kraken2_concordance"]["spearman_bt2rpm_vs_k2rpm"],confirmed_fraction=au["kraken2_concordance"]["overall_concordance"],breadth=au["breadth_26695"]),
    }


def reports(meta):
    rows=sorted(ROWS,key=lambda r:(r["check"],r["analysis"]))
    ROWS[:]=rows
    date_counts={d:sum(r['computed_date']==d for r in rows) for d in sorted({r['computed_date'] for r in rows})}
    timing="; ".join(f"{d}: {n} intervals" for d,n in date_counts.items())
    report=["# Topic B post-hoc CI computation", "", LABEL+".", "",
        "Actual computation dates (Asia/Shanghai): "+timing+". The first complete output was written at "+meta.get("computation_completed_at","the recorded run completion time")+". The approved 12 Nov analysis date is retained separately from the actual computation dates; the E2a calculation is not backdated.", "",
        "Part 1 only. These intervals were not pre-specified and do not change the original analysis tiers, P values, q values or multiplicity families. No manuscript Results/Methods revision is written.", "",
        "## Method and reproduction", "",
        "Two-sided 95% percentile intervals use 10,000 resamples and seed 20261112. Each analysis restarts the NumPy generator at that seed; the R E2 analysis uses R's seed 20261112 with fixed per-resample sub-seeds. Pairs are resampled as patients, including zero differences. The overall A1 RNA-16S correlation keeps all observed rows of each patient together. E12 resamples subjects with their group and outcome together (ordinary subject bootstrap), recalculating Jmax in each resample. NumPy percentile/R quantile use linear interpolation (R type 7). All selected point estimates must reproduce before their bootstrap is accepted.", "",
        "Saved score/deconvolution associations are conditional on the original constructed scores and BayesPrism expression estimates. BayesPrism is not refit. E1 repeats the original prevalence/rCLR preprocessing; E2 repeats score construction, expression filtering, TMM, voom, limma and the installed fgsea NES branch in each patient resample. Its vectorized TMM/voom/weighted-fit algebra is validated against the original R path on the full data and five deterministic bootstrap draws (t statistics and NES tolerance 1e-8); installed R lowess, limma span/variance moderation and fgsea code are reused. Original tests and P/q values are copied unchanged; the bootstrap is for estimation only. R2 is nonnegative and JT's null is J/Jmax=0.5, so excluding zero is not an appropriate significance check for these effects.", "",
        "## CHECK decisions", "", "| CHECK | Location | Decision |", "|---|---|---|"]
    report += [f"| {a} | {b} | {c} |" for a,b,c in DECISIONS]
    support=meta["supporting_effects_without_CI"]
    report += ["", "Supporting effects retained without new CIs (copied, not refit):", "",
        f"- C04 F2 Gemella: rho={fmt(support['C04']['reported_Gemella']['rho'])}, P={fmt(support['C04']['reported_Gemella']['P'])}.",
        f"- C07 E2b paired change: rho={fmt(support['C07']['paired']['rho'])}, P={fmt(support['C07']['paired']['P'])}; tumour: rho={fmt(support['C07']['tumour']['rho'])}, P={fmt(support['C07']['tumour']['P'])}.",
        f"- C08 unplotted tumour-only MUC5AC: rho={fmt(support['C08_tumour_context']['tumour_MUC5AC']['rho'])}, P={fmt(support['C08_tumour_context']['tumour_MUC5AC']['P'])}.",
        f"- C09 normal HP versus delta M1: rho={fmt(support['C09']['normal_HP_delta_M1']['rho'])}, P={fmt(support['C09']['normal_HP_delta_M1']['P'])}; normal HP versus IM: rho={fmt(support['C09']['normal_HP_IM']['rho'])}, P={fmt(support['C09']['normal_HP_IM']['P'])}.",
        f"- C12 classifier rho={fmt(support['C12']['classifier_rho'])}, confirmed fraction={fmt(support['C12']['confirmed_fraction'])}; classifier P remains unavailable. Breadth comparison P={fmt(support['C12']['breadth']['MW_breadthratio_N_vs_T_P'])}; original tissue medians are retained in JSON, with no invented effect contrast."]
    report += ["", "## Point reproduction and all new intervals", "", "| CHECK | Analysis | n | Original | Reproduced | Post-hoc 95% CI | Reproduction |", "|---|---|---:|---:|---:|---|---|"]
    report += [f"| {r['check']} | {r['analysis']} | {r['n']} | {fmt(r['original_estimate'])} | {fmt(r['estimate'])} | [{fmt(r['ci_low'])}, {fmt(r['ci_high'])}] | PASS (absolute error {r['absolute_reproduction_error']:.2g}) |" for r in rows]
    report += ["", "C02 is in proportions (multiply by 100 for percentage points). C15 is median paired T-N log2(epithelial CPM+1). C18 JT is J/Jmax; raw J intervals are also stored in JSON/CSV. All other correlations are Spearman rho; C11 is paired rank-biserial r, C03 is R2 and C06 is NES.", "", "## C01: Korean 16S patient accounting", ""]
    a=meta["C01"]
    report += [a["conclusion"],"",
        "The ENA manifest has 736 rows: 384 MiSeq WGS-labelled 16S runs plus 352 RNA-seq runs. The 384-run manifest and SILVA per-sample table have identical run sets. Removing the final two characters, as in `scripts/k16_analysis.py`, handles both `_N`/`_T` and `-N`/`-T` and yields 192 patients, each with exactly one N and one T. The RNA subset has 174 complete pairs and four tumour-only patients: "+", ".join(a["rna_tumour_only"])+". These are RNA missing normals, not 16S singletons.","",
        "Exact reconstruction of the erroneous 206: `library_name.str.replace(r'_[NT]$', '', regex=True)` collapses 356 underscore-suffix runs to 178 shared patients but leaves the 28 hyphen-suffix library IDs for "+", ".join(a["sixteen_s_only_patients"])+" distinct (for example GC13-N and GC13-T). Thus 178+28=206; correctly collapsing the latter gives 178+14=192. An underscore-only parser explains the exact discrepancy. The historical producer is not recoverable from retained scripts/logs, so this exact reconstruction is not represented as recovered source code. The cited sensitivity script does not write dedup_failbias.json. The verified counts resolve the manuscript denominator; the erroneous original field is left untouched.","",
        "E1's separate >=5000-read rule excludes four otherwise complete SILVA pairs ("+", ".join(a["E1_depth_excluded_pairs"])+"), leaving 188. This depth subset does not explain the erroneous total of 206.","", "## C14: BayesPrism deterministic exclusions", "",
        "Read-only rerun of the original mixture cleanup and protein-coding selection, followed by installed BayesPrism 2.2.3 `filter.bulk.outlier(outlier.cut=0.01, outlier.fraction=0.1)`. A gene is removed when its fraction of a sample's post-cleanup counts exceeds 1% in more than 10% of samples. Neither `new.prism` nor `run.prism` is called. Final shared gene sets are checked against each saved `bp_*.rds` prism mixture.",""]
    for cohort in ["korea","tcga"]:
        x=meta["C14"][cohort]
        report += [f"**{cohort}:** {x['n_outlier_genes']} outliers: "+", ".join(x["outlier_genes"])+".","",
            f"{x['n_samples']} samples; {x['n_after_protein_coding']} genes before outlier removal; {x['n_shared_retained']} retained reference-mixture genes. Stored mixture gene-set comparison: PASS.",""]
    report += ["MUC5AC and MUC6 are explicitly in the Korean outlier list. Gastric differentiation retains 5/10 genes: PGC, PGA3, LIPF and MUC6 were removed as outliers; CHIA was absent from the marker set. M3 retains 5/8: MUC5AC and MUC6 were outliers; B3GNT6 was absent from the marker set. Thus the failed M1/M3 epithelial scores have specific, recoverable causes. ST3GAL4 is also absent from the marker set. TCGA KRT4 is one of eight bulk outliers but was not in the marker reference, so only seven of those eight remove shared marker genes.","",
        "## CIs including the null and interpretation",""]
    affirmative=[r for r in rows if r["affirmative_claim"] and r["includes_null"]]
    report += ["Affirmative claims whose new interval includes the null: "+("; ".join(f"{r['analysis']} [{fmt(r['ci_low'])}, {fmt(r['ci_high'])}]" for r in affirmative) if affirmative else "none")+".",""]
    report += ["For E2a, this includes effects selected as significant and displayed in main Figure 5E. The five named main-text E2a estimates (MYC_TARGETS_V1 and the four frozen sets) all retain intervals excluding zero; four other plotted Hallmark sets include zero. Keep their exploratory status and original P/q values visible.",""]
    report += ["All intervals including their estimand's null: "+("; ".join(r["analysis"] for r in rows if r["includes_null"]))+".","",
        "The PRJNA310127 no-HP R2 point estimate (0.00843350) lies below its percentile interval [0.0102321, 0.0223633]. This is a limitation of the raw percentile interval for this bounded statistic, not a replacement for its original permutation P=0.0714. No supported tissue-location association is newly asserted.","",
        "A CI and an original P value can differ near a threshold because the bootstrap and original test use different sampling approximations. No new significance declaration or changed gene verdict is made in Part 1.","",
        "## Provenance and limits", "",
        "Full-precision estimates, reproduction errors, original P/q, exact source paths, resampling units, valid/invalid replicate counts and input SHA-256 hashes are in posthoc_ci.json; all interval rows are also in posthoc_ci.csv. Only the five authorised new files are written. Existing script/result/draft file sizes and nanosecond mtimes are checked before and after computation. The original C01 producer remains unavailable; this is the only historical provenance limitation in that accounting finding.", ""]
    summary="\n".join(report)
    (OUT/"POSTHOC_CI_SUMMARY.md").write_text(summary)
    note=summary.replace("# Topic B post-hoc CI computation","# Results skeleton v1b resolution notes — computation only",1)
    (ROOT/"manuscript/drafts/skeleton_v1b_resolution.md").write_text(note)
    csv_rows=[]
    for r in rows:
        csv_rows.append({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in r.items()})
    pd.DataFrame(csv_rows).to_csv(OUT/"posthoc_ci.csv",index=False)


def main():
    baseline=snapshot()
    for f in ["manuscript/drafts/results_skeleton_v1.md","manuscript/figure_allocation.md","results/figures/assembled/README.md"]:path(f)
    meta=dict(status="in_progress",label=LABEL,seed=SEED,n_resamples=NBOOT,method="percentile bootstrap",
              approved_analysis_date="2026-11-12",computation_started_at=datetime.now(ZONE).isoformat(),
              decisions=[dict(check=a,location=b,decision=c) for a,b,c in DECISIONS])
    if "--resume" in sys.argv and (OUT/"posthoc_ci.json").exists():
        saved=json.loads((OUT/"posthoc_ci.json").read_text())
        if saved.get("status")=="complete":
            for name,digest in saved["input_sha256"].items():
                assert hashlib.file_digest(path(name).open('rb'),'sha256').hexdigest()==digest, "Cached input changed: "+name
            ROWS.extend(saved["intervals"])
            meta={key:value for key,value in saved.items() if key!="intervals"}
            meta["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
            reports(meta);write_checkpoint(meta)
            assert snapshot()==baseline
            progress("Completed outputs verified and regenerated without repeating the bootstrap.")
            return
        ROWS.extend(r for r in saved["intervals"] if r["check"] not in ["C17","C18"])
        if "C14" in saved:meta["C14"]=saved["C14"]
        if "E2_bootstrap_checkpoint" in saved:meta["E2_bootstrap_checkpoint"]=saved["E2_bootstrap_checkpoint"]
        for r in ROWS:
            for f in [r["source"].split("::")[0],r["original_code"],*r["inputs"]]:path(f)
    meta["C01"]=accounting();progress("C01 verified source counts; historical producer not retained")
    k=read_csv("results/core/korea_rna_per_sample.csv",dtype={"pid":str})
    for check,fn in [("C02",lambda:a1(k)),("C03",community),("C08",lambda:e3_hm(k)),("C11",audit),("C15",lambda:e9(k)),("C17",a4),("C18",e12)]:
        if not any(r["check"]==check for r in ROWS):fn();write_checkpoint(meta)
    for f in ["scripts/s10_e7_bayesprism.R","results/step10/E7/marker_genes.txt","results/step10/E7/input/korea_bulk_counts.csv.gz",
              "data/host/star/stad_star_se.rds","results/core/hm_deltas_tcga.csv","data/external/spectrum2023/stad_helicobacter_per_sample.csv",
              "results/step10/E7/bp_korea.rds","results/step10/E7/bp_tcga.rds","results/step10/E7/E7_results.json","results/step10/E7/log/run_nc4.log"]:path(f)
    if "C14" not in meta:meta["C14"]=run_r(BAYES_R);write_checkpoint(meta)
    meta["E2_reproduction"]=e2(meta)
    meta["supporting_effects_without_CI"]=supporting_results()
    meta["input_sha256"]={str(p.relative_to(ROOT) if p.is_relative_to(ROOT) else p):hashlib.file_digest(p.open('rb'),'sha256').hexdigest() for p in sorted(INPUTS)}
    after=snapshot()
    assert after==baseline,"An existing project file changed or an unauthorised file was created"
    meta["file_integrity"]={"unchanged_existing_file_count":len(baseline),"changed_existing_files":[],"unexpected_new_files":[]}
    meta["status"]="complete";meta["n_intervals"]=len(ROWS)
    meta["computation_completed_at"]=datetime.now(ZONE).isoformat()
    meta["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    reports(meta);write_checkpoint(meta)
    assert snapshot()==baseline
    progress(f"Completed {len(ROWS)} intervals; all point estimates reproduced; existing files unchanged.")


if __name__=="__main__":
    main()
