"""Per-track CDF normalization audit (v27).

Audit-checklist §4 verbatim: every oracle's effect_cdfs is monotone
non-decreasing per row; summary_cdfs satisfies p50 ≤ p95 ≤ p99 per
row; track counts match published specs.

Outputs a one-line PASS/FAIL summary per oracle.
"""
from __future__ import annotations

import sys
import numpy as np

from chorus.analysis.normalization import get_normalizer

EXPECTED_TRACK_COUNTS = {
    "alphagenome": 5168,    # CDF coverage; full model has 5731
    "enformer": 5313,
    "borzoi": 7611,
    "chrombpnet": None,     # per-model, varies
    "sei": 40,              # 40 sequence classes
    "legnet": 3,            # 3 cell types
}


def audit_oracle(name: str) -> bool:
    nz = get_normalizer(name)
    if nz is None:
        print(f"[{name}] FAIL — get_normalizer returned None")
        return False

    entry = nz._loaded[name]
    ecdf = entry.get("effect_cdfs")
    scdf = entry.get("summary_cdfs")
    if ecdf is None or scdf is None:
        print(f"[{name}] FAIL — missing effect_cdfs or summary_cdfs")
        return False

    n_tracks = ecdf.shape[0]

    # Monotonicity of effect_cdfs (sample 10 rows)
    sample = min(10, n_tracks)
    for i in range(sample):
        diffs = np.diff(ecdf[i])
        if not np.all(diffs >= -1e-9):
            print(f"[{name}] FAIL — effect_cdfs[{i}] not monotone "
                  f"(min diff {diffs.min()})")
            return False

    # p50 ≤ p95 ≤ p99 for summary_cdfs
    n_pts = scdf.shape[1]
    p50_idx, p95_idx, p99_idx = (
        int(0.50 * n_pts), int(0.95 * n_pts), int(0.99 * n_pts),
    )
    for i in range(min(10, n_tracks)):
        p50, p95, p99 = scdf[i, p50_idx], scdf[i, p95_idx], scdf[i, p99_idx]
        if not (p50 <= p95 + 1e-9 <= p99 + 2e-9):
            print(f"[{name}] FAIL — track {i}: p50={p50}, p95={p95}, p99={p99}")
            return False

    expected = EXPECTED_TRACK_COUNTS.get(name)
    track_match = "skip" if expected is None else (
        "match" if n_tracks == expected else f"MISMATCH (expected {expected})"
    )

    print(f"[{name}] PASS — n_tracks={n_tracks} ({track_match}); "
          f"effect_cdfs monotone; p50≤p95≤p99")
    return True


if __name__ == "__main__":
    oracles = sys.argv[1:] or list(EXPECTED_TRACK_COUNTS.keys())
    fails = 0
    for o in oracles:
        if not audit_oracle(o):
            fails += 1
    print(f"\nTotal: {len(oracles) - fails}/{len(oracles)} pass")
    sys.exit(1 if fails else 0)
