# Build GDC open-access manifests for TCGA-STAD masked MAF and 450K beta values
import requests, json, pandas as pd
API="https://api.gdc.cancer.gov/files"
def q(dt, extra=[]):
    filt={"op":"and","content":[
      {"op":"in","content":{"field":"cases.project.project_id","value":["TCGA-STAD"]}},
      {"op":"in","content":{"field":"data_type","value":[dt]}},
      {"op":"in","content":{"field":"access","value":["open"]}}]+extra}
    fields="file_id,file_name,file_size,md5sum,platform,analysis.workflow_type,cases.samples.submitter_id,cases.samples.sample_type,cases.samples.portions.analytes.aliquots.submitter_id,cases.submitter_id"
    r=requests.post(API,json={"filters":filt,"fields":fields,"size":5000,"format":"JSON"},timeout=120).json()
    rows=[]
    for h in r["data"]["hits"]:
        for c in h["cases"]:
            for s in c["samples"]:
                al=[a["submitter_id"] for p in s.get("portions",[]) for an in p.get("analytes",[]) for a in an.get("aliquots",[])]
                rows.append(dict(file_id=h["file_id"],file_name=h["file_name"],size=h["file_size"],md5=h["md5sum"],platform=h.get("platform"),
                  workflow=h.get("analysis",{}).get("workflow_type"),case=c["submitter_id"],sample=s["submitter_id"],sample_type=s["sample_type"],aliquots=";".join(al)))
    return pd.DataFrame(rows)
m=q("Masked Somatic Mutation"); m.to_csv("data/a6/manifest_maf.csv",index=False)
b=q("Methylation Beta Value",[{"op":"in","content":{"field":"platform","value":["Illumina Human Methylation 450"]}}]); b.to_csv("data/a6/manifest_meth450.csv",index=False)
for n,d in [("maf",m),("meth",b)]:
    print(n,d.file_id.nunique(),"files",d.drop_duplicates("file_id")["size"].sum()/1e9,"GB");print(d.sample_type.value_counts().to_dict());print(d.workflow.value_counts().to_dict())
