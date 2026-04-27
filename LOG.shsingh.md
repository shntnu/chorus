# LOG.shsingh.md

Day-by-day exploration notes for this fork. Oldest first.
Standing context lives in `CLAUDE.shsingh.md`. Promote actionable
findings to `../pign-cdg/docs/log.md` with a link back here.

## 2026-04-27 - Fork setup, Enformer smoke test on Mac + Oppy

- Forked `pinellolab/chorus`. Added a Nix flake (micromamba-based,
  since chorus's per-oracle envs are conda-managed and TF/PyTorch/JAX
  don't coexist in one env). `.envrc` with `use flake`. `Justfile`
  wraps `chorus setup` and a `just demo` smoke test.
- macOS: env created and Enformer downloaded fine; cancelled the
  CDF-backgrounds download mid-stream (planning to move to Oppy
  anyway).
- Oppy (NixOS, 4x H100 NVL): full `just setup` ran clean. Enformer
  weights, hg38, and CDF backgrounds all pulled. `just demo` returned
  real numbers (WT DNase signal at chr11:5,247,000-5,248,000 = 0.468;
  3 alt alleles scored at chr11:5,247,500). Linux CUDA stack
  (cuDNN 8, cuBLAS, etc.) installed via the per-oracle env's pip
  deps - so the env was platform-aware (TF metal on Mac, CUDA on
  Linux).
- Did not install AlphaGenome yet. Next time: add it (need
  `HF_TOKEN` for the gated model), then run a real PIGN-shaped query.

## 2026-04-27 - AlphaGenome install on Oppy

- Auth: `micromamba run -n chorus hf auth login`. Note
  `huggingface-cli` is deprecated, replaced by `hf` (chorus's
  `_tokens.py` calls `huggingface_hub.login()` directly so the binary
  rename doesn't matter for chorus itself, but the README's hint
  command is stale). Token persisted to
  `~/.cache/huggingface/token`.
- First `just oracle alphagenome` failed at the weight-prefetch step:
  403 on `google/alphagenome-all-folds`. Token was fine
  (fineGrained, scope: "Read access to public gated repos you can
  access") - my HF account just hadn't accepted Google's gate. Fixed
  by clicking through at
  <https://huggingface.co/google/alphagenome-all-folds>; granted
  immediately. Re-running was idempotent on env (already built,
  6.3 GB) and pulled weights in ~30 s
  (~700 MB at `~/.cache/huggingface/hub/models--google--alphagenome-all-folds`).
- chorus rewrites `jax[cpu]` → `jax[cuda12]` based on detected
  hardware during pip install (visible in install log line
  "Installing pip packages: jax[cuda12], ..."). The yml's `jax[cpu]`
  is just a platform-neutral starting point. Earlier worry that
  AlphaGenome would be CPU-only on Oppy was unfounded.
- Smoke test, parity with the Enformer `just demo` (β-globin locus,
  K562 DNase, chr11:5,247,000-5,248,000):
  - Track id: `DNASE/EFO:0002067 DNase-seq/.` (only 1 K562 DNase
    track in AlphaGenome's 5,731-track set; cell-line-level assays
    use EFO ontology, primary-cell-type assays use CL).
  - Model load: 31 s. WT prediction: 57 s. Variant scan
    (ref C → A/G/T at chr11:5,247,500): 232 s for 4 alleles
    (~58 s each).
  - WT mean DNase: 0.045. Per-allele deltas: <1e-4. Both numbers
    are misleading on their own — AlphaGenome forces 1 MB context,
    so the mean is over 1 MB (vs. Enformer's 1 kb-window mean of
    0.468 — not comparable). Max signal in window: ~50, consistent
    with strong DHSs in the β-globin LCR. Variant deltas are tiny
    because a single base in 1 MB barely moves the global mean —
    real variant scoring needs a local-windowed strategy, not the
    global mean.
- GPU utilization unverified. chorus logs only said
  "Device: auto-detect"; I never watched `nvidia-smi` during the
  run. Confirm next time with `nvidia-smi dmon`.

## Open threads

- First real query: PIGN expression across chorus's CAGE / RNA-seq
  tracks - which cell types light up, how does that line up with
  what's known about MCAHS1 affected tissues?
- Pull regulatory-region PIGN variants from ClinVar (most are coding,
  but check) and score the small set if any exist.
- Variant scoring strategy: chorus's per-track `preferred_scoring_strategy`
  defaults to `"mean"` for AlphaGenome, but for variant effects in a
  1 MB context that's too dilute. Check what local-window aggregation
  chorus exposes (likely on `OraclePrediction` or a helper) before
  scoring real PIGN variants.
- `nvidia-smi dmon` during AlphaGenome inference to confirm GPU
  utilization (currently inferred from timing only).
