import subprocess, time, json, urllib.request, urllib.error, sys
P=8099; B=f"http://127.0.0.1:{P}"
srv=subprocess.Popen(["./crispr","serve",str(P)],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
time.sleep(0.5)
fails=0
def call(path,data=None):
    req=urllib.request.Request(B+path,data=data,method="POST" if data is not None else "GET")
    try:
        with urllib.request.urlopen(req,timeout=120) as r: return r.status,r.read()
    except urllib.error.HTTPError as e: return e.code,e.read()
def js(path,data=None,expect=200):
    global fails
    st,b=call(path,data)
    try: j=json.loads(b)
    except Exception as e:
        print("  FAIL not JSON:",path,b[:80]); fails+=1; return None
    ok=(st==expect)
    if not ok: fails+=1
    print(("  ok  " if ok else "  FAIL"),st,path[:70])
    return j
def chk(c,m):
    global fails
    print("   ",("pass" if c else "FAIL"),m)
    if not c: fails+=1
try:
    j=js("/api/status"); chk(j["loaded"] is False,"status before load: loaded=false")
    j=js("/api/find",expect=400); chk("error" in j,"find without genome gives plain error: "+j["error"][:60])
    j=js("/api/synth?n=5",expect=400); chk("error" in j,"synth n=5 rejected")
    j=js("/api/synth?n=300000"); chk(j["loaded"] and j["n"]==300000 and j["planted"]["pos"]==100000,"synth 300k loaded, planted at 100000")
    j=js("/api/find"); chk(j["agree"] and j["pattern16"]["agree"] and len(j["matchers"])==3,"find: 3 matchers agree (%d targets)"%j["matchers"][0]["targets"])
    j=js("/api/rank"); chk(j["ntop"] if False else len(j["top"])==10 and any(w["pos"]==100000 for w in j["worst"]),"rank default region: top 10, planted among worst")
    j=js("/api/rank?from=0&to=40000",expect=400); chk("30,000" in j["error"],"rank region >30000 rejected: "+j["error"][:50])
    j=js("/api/rank?from=100&to=90",expect=400); chk("error" in j,"rank to<from rejected")
    j=js("/api/rank?K=abc",expect=400); chk("error" in j,"rank bad K rejected")
    j=js("/api/rank?from=90000&to=110000&K=5&k=5"); chk(len(j["top"])==5 and j["k"]==5,"rank k=5 (index rebuilt: %.1f ms)"%j["index_ms"])
    j=js("/api/hits?guide=GAGTCCGAGCAGAAGAAGAA&k=3"); chk(j["count"]==13 and sum(1 for h in j["hits"] if h["mm"]==0)==1,"hits: 13 sites (1 intended + 12 decoys)")
    chk(all(len(h["site"])==23 and len(h["mask"])==20 and h["mask"].count("1")==h["mm"] for h in j["hits"]),"hits: site=23 letters, mask matches mm count")
    j=js("/api/hits?guide=ACGT",expect=400); chk("20 letters" in j["error"],"hits bad guide rejected")
    j=js("/api/offtarget?k=3&n=10"); chk(j["identical"]==j["of"]==10,"offtarget: identical %d of %d, speedup %.0fx per query, %.1fx with build"%(j["identical"],j["of"],j["speedup_query"],j["speedup_with_build"]))
    j=js("/api/offtarget?k=9",expect=400); chk("error" in j,"offtarget k=9 rejected")
    j=js("/api/cpg"); chk(j["count"]==3 and len(j["bins"])==200,"cpg: 3 islands, 200 GC bins")
    j=js("/api/align?len=1000"); chk(j["equal"] is True and j["hirschberg"]["bytes"]<j["nw"]["bytes"]/10,"align len=1000: scores equal (%d), memory %d vs %d bytes"%(j["nw"]["score"],j["nw"]["bytes"],j["hirschberg"]["bytes"]))
    j=js("/api/align?a=GATTACA&b=GCATGCU",expect=400); chk("error" in j,"align bad letters rejected")
    j=js("/api/align?a=GATTACA&b=GCATGCA"); chk(j["equal"] is True and j["m"]==7,"align custom input")
    j=js("/api/align?a=GATTACA",expect=400); chk("error" in j,"align with only one sequence rejected")
    j=js("/api/align?len=8000"); chk(j["equal"] is True,"align len=8000: equal (%d), NW %d bytes vs Hirschberg %d"%(j["nw"]["score"],j["nw"]["bytes"],j["hirschberg"]["bytes"]))
    j=js("/api/vault?bits=64"); chk(all(c["pass"] for c in j["checks"]) and j["restored"]==j["original"],"vault: all 3 checks pass (%d -> %d -> %d bytes)"%(j["original"],j["compressed"],j["encrypted"]))
    j=js("/api/vault?bits=10",expect=400); chk("error" in j,"vault bits=10 rejected")
    j=js("/api/attack"); chk(all(r["recovered"]==r["trials"] for r in j["rows"]) and len(j["extrapolation"])==4,"attack: %d rows, all keys recovered; 1024-bit: ~10^%.0f years"%(len(j["rows"]),j["extrapolation"][3]["log10_years"]))
    j=js("/api/bench"); chk(len(j["genomes"])==4 and len(j["align"])==4,"bench: 4 genome sizes, 4 alignment sizes")
    # uploads
    import os
    subprocess.run(["./crispr","gen","60000","/tmp/up.fa","5"],capture_output=True)
    fa=open("/tmp/up.fa","rb").read()
    j=js("/api/upload?name=my%20genome.fa",fa); chk(j["loaded"] and j["n"]==60000 and j["name"]=="my genome.fa" and j["planted"] is None,"upload FASTA: 60000 bases, name kept, no planted info")
    j=js("/api/upload?name=junk.txt",b"this is not a genome\nat all\n",expect=400); chk("FASTA" in j["error"],"upload garbage: "+j["error"][:70])
    j=js("/api/upload",b"",expect=400); chk("error" in j,"upload empty body: "+j["error"])
    j=js("/api/upload",b">x\r\nNNNNNNNN\r\n",expect=400); chk("A, C, G or T" in j["error"],"upload FASTA with no ACGT: "+j["error"])
    j=js("/api/status"); chk(j["name"]=="my genome.fa","genome still the last good upload after failed uploads")
    st,b=call("/api/upload"); chk(st==400,"GET on upload -> 400")
    st,b=call("/nope"); chk(st==404,"unknown path -> 404")
    st,b=call("/"); chk(st==200 and b"<" in b,"/ serves web/index.html")
    # big upload (4.6 MB)
    subprocess.run(["./crispr","gen","4600000","/tmp/big.fa","3"],capture_output=True)
    j=js("/api/upload?name=big.fa",open("/tmp/big.fa","rb").read()); chk(j["n"]==4600000,"upload 4.6 MB FASTA")
finally:
    srv.terminate(); srv.wait(timeout=5)
print("\nAPI test:", "ALL PASSED" if fails==0 else f"{fails} FAILED")
sys.exit(fails!=0)
