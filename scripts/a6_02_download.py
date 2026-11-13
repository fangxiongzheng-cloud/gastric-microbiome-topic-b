# Resumable GDC download (open access) with md5 check. Usage: python a6_02_download.py maf|meth
import sys, os, hashlib, requests, pandas as pd, time
kind=sys.argv[1]
hp=pd.read_csv("results/stage1/hp_calls.csv")
hp=hp[(hp.tissue=="T")&hp.in380&hp.hp_consensus.isin(["pos","neg"])]
cases=set(hp.case)
m=pd.read_csv(f"data/a6/manifest_{'maf' if kind=='maf' else 'meth450'}.csv")
if kind=="meth": m=m[(m.sample_type=="Primary Tumor")&m.case.isin(cases)]
m=m.drop_duplicates("file_id")
out=f"data/a6/{kind}"; os.makedirs(out,exist_ok=True)
def md5(p):
    h=hashlib.md5(); f=open(p,"rb")
    for ch in iter(lambda:f.read(1<<20),b""): h.update(ch)
    return h.hexdigest()
ok=0
for i,r in enumerate(m.itertuples()):
    p=f"{out}/{r.file_id}__{r.file_name}"
    if os.path.exists(p) and md5(p)==r.md5: ok+=1; continue
    for t in range(5):
        try:
            with requests.get(f"https://api.gdc.cancer.gov/data/{r.file_id}",stream=True,timeout=300) as g:
                g.raise_for_status()
                with open(p+".part","wb") as f:
                    for ch in g.iter_content(1<<20): f.write(ch)
            if md5(p+".part")==r.md5: os.replace(p+".part",p); ok+=1; break
        except Exception as e: print("retry",r.file_id,e,flush=True); time.sleep(10)
    if i%25==0: print(kind,i,len(m),flush=True)
print("DONE",kind,ok,"/",len(m),flush=True)
