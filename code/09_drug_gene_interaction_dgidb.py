# -*- coding: utf-8 -*-
"""
用真实的 DGIdb GraphQL API 查询 7 个标志物的药物-基因相互作用。

DGIdb 已从 REST 迁移到 GraphQL（dgidb.org/api/v2 不再有效），因此本脚本使用
GraphQL 端点查询，并把完整原始结果与按规则过滤后的结果分别落盘，保证可复现、
可核查。

输出：
  results/dgidb_raw_interactions.csv        全部原始记录
  results/dgidb_drug_gene_interactions.csv  过滤后（获批药物 且 interactionScore >= 1）
"""

import csv
import json
import os
import urllib.request

PROXY = "http://127.0.0.1:7897"
GENES = ["NCKAP1L", "THSD7A", "MYH10", "PDLIM1", "IL1B", "FLNB", "SLC3A2"]
URL = "https://dgidb.org/api/graphql"
ROOT = r"C:\Users\weihong\Desktop\生物竞赛"
OUTDIR = os.path.join(ROOT, "results")
SCORE_CUTOFF = 1.0
SCORE_FLOOR = 2.0
TYPED_KEEP = ("inhibitor", "antibody", "immunotherapy", "blocker", "negative modulator")

QUERY = """
query ($names: [String!]!) {
  genes(names: $names) {
    nodes {
      name
      interactions {
        drug { name approved }
        interactionScore
        interactionTypes { type directionality }
        sources { sourceDbName }
      }
    }
  }
}
"""


def fetch():
    body = json.dumps({"query": QUERY, "variables": {"names": GENES}}).encode("utf-8")
    req = urllib.request.Request(
        URL, data=body,
        headers={"Content-Type": "application/json", "User-Agent": "codex-dgidb-client"},
    )
    opener = urllib.request.build_opener(
        urllib.request.ProxyHandler({"http": PROXY, "https": PROXY})
    )
    with opener.open(req, timeout=180) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    if "errors" in payload:
        raise RuntimeError(payload["errors"])
    return payload["data"]["genes"]["nodes"]


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    nodes = fetch()

    rows = []
    for node in nodes:
        gene = node["name"]
        for ix in node["interactions"]:
            types = ix.get("interactionTypes") or []
            type_str = ";".join(
                "%s%s" % (t["type"], ("(%s)" % t["directionality"]) if t.get("directionality") else "")
                for t in types
            )
            sources = [s["sourceDbName"] for s in (ix.get("sources") or [])]
            rows.append({
                "Target_Gene": gene,
                "Drug": ix["drug"]["name"],
                "Approved": ix["drug"]["approved"],
                "InteractionScore": round(float(ix["interactionScore"]), 4),
                "InteractionType": type_str,
                "Sources": ";".join(sources),
                "N_Sources": len(sources),
            })

    rows.sort(key=lambda r: (r["Target_Gene"], -r["InteractionScore"]))

    raw_path = os.path.join(OUTDIR, "dgidb_raw_interactions.csv")
    with open(raw_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # 过滤规则（论文中需原文说明）：
    #   A. 已获批药物 且 interactionScore >= 2.0；或
    #   B. 带明确实验证据的相互作用类型（抑制剂/抗体/免疫治疗等）且被 >= 2 个独立来源数据库支持
    kept, seen = [], set()
    for r in rows:
        typed = any(t in r["InteractionType"].lower() for t in TYPED_KEEP)
        clause_a = r["Approved"] and r["InteractionScore"] >= SCORE_FLOOR
        clause_b = typed and r["N_Sources"] >= 2
        if not (clause_a or clause_b):
            continue
        key = (r["Target_Gene"], r["Drug"])
        if key in seen:
            continue
        seen.add(key)
        kept.append(r)
    kept.sort(key=lambda r: (r["Target_Gene"], -r["InteractionScore"]))

    filt_path = os.path.join(OUTDIR, "dgidb_drug_gene_interactions.csv")
    with open(filt_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(kept)

    print("原始记录: %d 条 -> %s" % (len(rows), raw_path))
    print("过滤后  : %d 条 -> %s" % (len(kept), filt_path))
    print()
    print("%-10s %8s %8s %8s" % ("gene", "raw", "kept", "maxScore"))
    for g in GENES:
        sub = [r for r in rows if r["Target_Gene"] == g]
        kk = [r for r in kept if r["Target_Gene"] == g]
        mx = max([r["InteractionScore"] for r in sub], default=0)
        print("%-10s %8d %8d %8.2f" % (g, len(sub), len(kk), mx))
    print()
    print("过滤后药物（按评分降序前 30）:")
    for r in kept[:30]:
        print("  %-10s %-32s %5.2f  src=%s" % (r["Target_Gene"], r["Drug"], r["InteractionScore"], r["N_Sources"]))


if __name__ == "__main__":
    main()
