"""
Follow-up (post hoc) analysis: does translation direction matter for the open-weight
group and for the closed-source group?  6 tests total (2 groups x 3 metrics).

Standalone: does NOT rewrite any of your existing tables. Writes one new file:
  tables/table_direction_by_group.csv
Requires: pandas, numpy, scipy
"""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

# ---- same config as your main script ----
CSV = "JOKES - Sheet5.csv"
METRICS = {"Final_Humour": "Humor", "Final_Fluency": "Fluency", "Final_Accuracy": "Accuracy"}
CLOSED = ["gpt-5.2-2025-12-11", "gemini2.5-flash", "claude-sonnet-4-5"]
OPEN = ["llama3.1:latest", "qwen3:8b", "deepseek-r1:latest"]
GROUPS = {"Open-weight": OPEN, "Closed-source": CLOSED}
OUT = Path("tables")
OUT.mkdir(exist_ok=True)


def bh(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    adj = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adj, 1)
    return out


def perm_indep_p(a, b):
    """Exact two-sided permutation test for a difference in means (independent groups)."""
    pooled = np.concatenate([a, b])
    n, N = len(a), len(pooled)
    obs = abs(a.mean() - b.mean())
    hits = tot = 0
    for c in combinations(range(N), n):
        mask = np.zeros(N, bool)
        mask[list(c)] = True
        hits += abs(pooled[mask].mean() - pooled[~mask].mean()) >= obs - 1e-12
        tot += 1
    return hits / tot


df = pd.read_csv(CSV)
df["Joke_UID"] = df["Direction"] + "_J" + df["Joke ID"].astype(str)
directions = list(df["Direction"].unique())
assert len(directions) == 2, f"Expected 2 directions, found {directions}"

rows = []
for gname, members in GROUPS.items():
    sub = df[df["Model"].isin(members)]
    for col, label in METRICS.items():
        # one value per joke = average of that group's 3 models on that joke
        per_joke = sub.groupby(["Direction", "Joke_UID"])[col].mean().reset_index()
        g1 = per_joke[per_joke["Direction"] == directions[0]][col].values
        g2 = per_joke[per_joke["Direction"] == directions[1]][col].values
        rows.append({
            "Group": gname, "Metric": label, "n1": len(g1), "n2": len(g2),
            f"Mean_{directions[0]}": g1.mean(), f"Mean_{directions[1]}": g2.mean(),
            "Diff_(dir0 - dir1)": g1.mean() - g2.mean(),
            "U": stats.mannwhitneyu(g1, g2, alternative="two-sided").statistic,
            "p_permutation": perm_indep_p(g1, g2),
        })
res = pd.DataFrame(rows)
res["p_BH"] = bh(res["p_permutation"].values)   # one family: these 6 tests
res.round(4).to_csv(OUT / "table_direction_by_group.csv", index=False)

print("=== DIRECTION EFFECT BY GROUP (post hoc; jokes = units; BH across 6 tests) ===")
print(res.round(4).to_string(index=False))
print(f"\nDone. Wrote {(OUT / 'table_direction_by_group.csv').resolve()}")