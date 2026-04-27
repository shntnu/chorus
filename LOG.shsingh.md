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

## Open threads

- AlphaGenome install on Oppy.
- First real query: PIGN expression across chorus's CAGE / RNA-seq
  tracks - which cell types light up, how does that line up with
  what's known about MCAHS1 affected tissues?
- Pull regulatory-region PIGN variants from ClinVar (most are coding,
  but check) and score the small set if any exist.
