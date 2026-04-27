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

## 2026-04-27 - GPU verification: was NOT engaged; libcuda.so path fix

- Did the verification check that the previous entry kicked down the
  road. Inside chorus-alphagenome: `jax.devices()` returned
  `[CpuDevice(id=0)]` with `cuInit(0) failed: error 303 ... unable to
  load CUDA libraries`. JAX silently fell back to CPU. So the earlier
  "AlphaGenome would be CPU-only ... was unfounded" sentence was
  itself wrong - it WAS CPU-only on Oppy, just for a different reason
  than I was worried about. `jax[cuda12]` was installed correctly;
  the runtime couldn't dlopen `libcuda.so`. "Device: auto-detect" in
  chorus logs is the chorus base process talking, not the per-oracle
  subprocess - unreliable.
- Cross-check on chorus-enformer (TF): same fault, different error
  ("Cannot dlopen some GPU libraries", plural). Systemic, not
  AlphaGenome-specific.
- Root cause: on NixOS, `libcuda.so` lives at
  `/run/opengl-driver/lib/libcuda.so` (symlink → nix-store
  `nvidia-dc-565.57.01` derivation). chorus's flake.nix `libList →
  LD_LIBRARY_PATH` did not include that path, so pip-installed CUDA
  wheels couldn't find it.
- Fix: one-line shellHook addition mirroring the validated pattern
  from `shntnu-neusis/templates/python-pixi/flake.nix:31`
  (Linux-guarded so the same flake works on Mac):
  `[ -d /run/opengl-driver/lib ] && export LD_LIBRARY_PATH=/run/opengl-driver/lib:$LD_LIBRARY_PATH`.
  Note from upstream of that template (MAINTENANCE_LOG line 760):
  "Validated pixi for GPU/RAPIDS workflows without FHS wrapper
  complexity" - the LD path is enough; no buildFHSEnv / nix-ld
  acrobatics needed for CUDA specifically. The path is set on the
  dev shell and inherited by anything spawned inside (micromamba,
  pixi, plain pip).
- Verified after fix, inside chorus-alphagenome env (with my
  LD_LIBRARY_PATH override matching the new flake):
  - `jax.devices()` → `[CudaDevice(id=0..3)]`, default backend=gpu.
  - 4Kx4K matmul: **66x speedup** (0.54 ms/call on GPU vs 35.56
    ms/call on CPU).
  - chorus's per-oracle subprocess inherits the path correctly
    (verified by running an introspection script through
    `oracle._env_runner.run_script_in_environment`); jax inside the
    subprocess also sees all 4 GPUs.
- Surprise: AlphaGenome's chorus-level prediction timing didn't
  improve (~64 s/call, same as CPU). Cause: `predict_template.py:48`
  calls `create_from_huggingface(fold, device=jax_device)` inside
  every subprocess, so each `oracle.predict(...)` does a full
  model-load + JAX/XLA compile + inference. The ~64 s is dominated
  by the load+compile (~50 s); actual GPU inference is the small
  remainder. Two back-to-back `oracle.predict()` calls both took
  64.4 s and 64.3 s - 1.0× speedup, not because GPU isn't working
  but because the compile cache doesn't survive subprocess death.
  This is a chorus architectural property (per-call subprocess
  isolation across oracle envs), not a fix bug. Real cohort work
  will need either chorus's batch interface (if any) or dropping
  down to AlphaGenome's API directly.
- Enformer (chorus-enformer): error after the fix changed from "DSO
  not loaded" to "Cannot dlopen some GPU libraries (plural)". TF
  2.13.1 predates the bundled-CUDA-wheels era (TF 2.15+ ships them);
  it needs cudatoolkit/cudnn/cublas installed externally. Not fixed
  by the LD_LIBRARY_PATH change. Out of scope for now.

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
- chorus per-prediction subprocess overhead (~60 s/call) hides GPU
  speedup for single calls because the model + XLA compile don't
  survive between subprocess invocations. For PIGN cohort scoring,
  check whether chorus has a batch / persistent-session mode, or if
  we'd need to drop down to alphagenome's API directly to amortize
  the compile.
- Enformer GPU: chorus-enformer pins TF 2.13.1, which predates the
  bundled-CUDA-wheels era. Either install cudatoolkit+cudnn into the
  env, or bump chorus-enformer.yml to TF 2.15+. Probably not worth
  chasing unless we actually need Enformer for the PIGN work -
  AlphaGenome supersedes it for our use case.
