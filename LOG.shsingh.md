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

## 2026-04-27 - Enformer GPU: thread closed, was a misdiagnosis

- Re-checked the parked "Enformer GPU broken" claim from the
  19:34 entry. It's wrong - the GPU path actually works.
- The CUDA-11 wheels are already in `environments/chorus-enformer.yml:42-49`
  (`nvidia-cublas-cu11`, `nvidia-cudnn-cu11>=8.6,<9.0`, `cuda_runtime`,
  `cufft`, `curand`, `cusolver`, `cusparse`). No external cudatoolkit
  install needed.
- `chorus.core.environment.runner._prepare_env` (runner.py:97-112)
  already stitches `nvidia/*/lib` onto `LD_LIBRARY_PATH` for any TF
  env. The comment in that code is exactly "TF <2.15 not auto-discovering
  nvidia-* pip packages".
- TF 2.13.1 inside chorus-enformer, when invoked through the
  EnvironmentRunner harness: `Built with CUDA: True`, four
  `PhysicalDevice('/physical_device:GPU:N')` visible.
- `just demo` runs end-to-end on GPU: WT mean DNase = 0.469 (matches
  the Mac CPU smoke 0.468 within rounding). nvidia-smi shows GPU 0
  jumping to 93 GB allocated (TF default greedy alloc). Wall clock
  is ~70 s end-to-end (re-verified 2026-04-27 evening): ~13 s for
  load + first predict, the rest dominated by `predict_variant_effect`
  spawning a fresh subprocess per alt allele. Earlier "~13 s" claim
  in this entry was timing only load + WT predict, not the full demo.
- Where the 19:34 misdiagnosis came from: same pattern as the
  AlphaGenome libcuda issue - the test was run *outside* the
  EnvironmentRunner harness (raw `mamba run -n chorus-enformer python
  -c ...`), which doesn't add the pip-installed nvidia libs to
  `LD_LIBRARY_PATH`. So both error messages ("DSO not loaded",
  "Cannot dlopen some GPU libraries") said "drivers missing" for
  the same wrong reason.
- One real but minor caveat: TF 2.13.1 logs `TensorFlow was not built
  with CUDA kernel binaries compatible with compute capability 9.0.
  CUDA kernels will be jit-compiled from PTX, which could take 30
  minutes or longer.` H100 is sm_90 (Hopper); TF 2.13's bundled
  kernels stop at sm_86. Practically, demo completed in seconds -
  PTX JIT happens but is fast for Enformer's kernel set. Worth
  remembering for fresh-process workloads, where the JIT cost is
  paid once per subprocess (and dies with it, same problem as
  AlphaGenome's compile cache).
- Lesson for future verification: any "GPU broken" claim in this
  repo needs to be reproduced *through* the EnvironmentRunner
  (`oracle.predict(...)` or a `just demo`-style harness call), not
  by a raw `mamba run`, because the runner does non-trivial env-prep
  that's invisible from a shell.

## 2026-04-27 — Why this fork exists, in plain language

Stepping back from the day's plumbing to write down, in one place, what
this exploration is for. Audience: someone in a PIGN-CDG family who
isn't a computational biologist and wants to know what I'm actually
doing on this side of the project.

### The bigger picture: pign-cdg

PIGN-CDG (also called MCAHS1) is a rare genetic condition. The PIGN
gene is one piece of an "anchor factory" inside every cell — a
~27-step assembly line that builds tiny molecular tethers (called GPI
anchors) which hold over 150 different proteins on the outside of the
cell, where they do their jobs. When PIGN doesn't work, those
proteins can't get tethered properly, and the downstream effects
touch many systems at once: the brain (seizures, developmental
delay), muscles (low tone), and others.

NS, my niece, was diagnosed with PIGN-CDG in 2024. She inherited
two different broken copies of PIGN — one from each parent. About
140 other families worldwide are known to be affected. The pign-cdg
project at the Broad Institute is the umbrella for the work my lab
and collaborators are doing to help; it has three main strands:

1. **Drug repurposing.** Perlara screened ~1,760 already-approved
   drugs and supplements against a "yeast avatar" of PIGN — yeast
   engineered to carry the human disease. Three rescued the yeast:
   ascorbyl palmitate (a vitamin C derivative), voriconazole (an
   antifungal), and ziprasidone (an antipsychotic). Some children,
   including NS for a stretch, have tried ascorbyl palmitate;
   results across families have been mixed.

2. **Variant characterization.** Different families carry different
   mutations in PIGN. We need to know, for each one: how broken is
   that specific version? My colleague RS has been analyzing images
   of cells in which 108 different PIGN variants — including NS's —
   have been overexpressed, characterizing how each variant alters
   cell appearance and protein localization. That work is the spine
   of the paper we're planning.

3. **This exploration: the DNA-regulation angle.** What this fork
   is about. See below.

The pign-cdg repo (`../pign-cdg`) is where the main analyses, the
patient brief, the collaborator map, and the paper plan all live.
That is the canonical place; this fork is a side-experiment.

### What I'm doing here, and what "chorus" / "LP" means

**LP** runs a computational biology lab at MGH / Harvard.
His group released **chorus** — an open-source software toolkit that
bundles together six different deep-learning models trained to
*read DNA*. Each model takes a stretch of DNA as input and predicts
what that DNA does in different cell types — for example, "this
piece of DNA acts as an on-switch for gene X in liver cells but is
silent in neurons." The six models in chorus (Enformer, AlphaGenome,
Borzoi, ChromBPNet, Sei, LegNet) were each built by different
research groups and trained on different data, so they're
complementary lenses on the same question: *what does this stretch
of DNA do?*

The first two strands of pign-cdg work are about the **protein**:
how broken is the PIGN protein when a particular mutation is in
place, and what drug might help. Chorus operates one level upstream
of that — given a stretch of DNA (anywhere in or around a gene), it
predicts molecular readouts that don't depend on knowing the protein
structure: how much of the gene is transcribed, where the chromatin
is open, where transcription factors bind, where splicing occurs.
Two specific questions it could help with:

- **Where is PIGN normally active?** Chorus can predict, across
  thousands of cell types, where PIGN is most expressed. If those
  predictions line up with the tissues most affected in MCAHS1
  (brain, muscle), it tells us we're using cell models that resemble
  the real disease tissue — useful for choosing the right cells in
  follow-up experiments.

- **Are there non-coding mutations we're missing?** Almost all known
  pathogenic PIGN mutations sit inside the protein-coding region.
  But for some of the ~140 families, no clear cause has been found
  in the coding region. Mutations in the regulatory DNA *around*
  PIGN — promoters, enhancers, splice sites — could explain some of
  those. Chorus is one of the only practical tools for scoring that
  class of mutation, and existing methods used elsewhere on this
  project (imaging, yeast functional assays, protein-folding
  predictors) don't cover it.

### Scope and limits

NS's two mutations are coding (a missense and a multi-exon deletion).
Chorus predicts gene-regulatory and splicing readouts from sequence —
not the protein-biochemistry effects of an amino-acid change (folding,
binding, enzyme activity). The central question for her variants —
does Arg95Gln disrupt PIGN's enzymatic role in GPI-anchor synthesis?
— sits in the protein-structural / variant-effect space (AlphaFold,
ESM, RS's image-based profiling). Chorus could still score these
variants for a secondary effect — e.g., does Arg95Gln happen to
disrupt a splice site or local expression — but that's a side check,
not the primary lens. What chorus covers that the other strands do
not:

- Tissue-of-interest selection: predicting where PIGN is most
  active across cell types, as a second lens alongside the
  imaging work.
- Non-coding mutations in the broader cohort: regulatory-region
  variants are not addressed by the imaging, yeast, or
  protein-folding tools currently in use.
- Cross-modality interpretation: combining imaging (RS),
  functional (yeast / VISTA), and regulatory (chorus) signals on
  a single variant.

Expected output: a supplementary section in the paper and any
leads from the cohort scan. The exploration takes a few weeks on
a server already in place.

### Where things stand today

Setup: the tool runs end-to-end on a 4-GPU server. NVIDIA driver
path issues on NixOS were resolved earlier today (see the prior
entries). The first scientific query — predicting PIGN's expression
landscape across cell types — is the next thing.

I have done **zero PIGN-specific analysis with chorus yet**. Today
was infrastructure. The biology starts tomorrow.

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
- chorus per-prediction subprocess overhead hides GPU speedup for
  single calls because the model + XLA/PTX-JIT compile don't survive
  between subprocess invocations. This is *all* oracles, not just
  AlphaGenome — confirmed by re-verification: `just demo` (Enformer,
  WT predict + 3-allele variant scan) takes ~70 s end-to-end on GPU
  because each alt allele re-pays the load cost. AlphaGenome is just
  the worst case (~60 s/call vs Enformer's ~13 s/call) because its
  model is bigger. For PIGN cohort scoring, check whether chorus has
  a batch / persistent-session mode, or if we'd need to drop down to
  the underlying oracle APIs directly to amortize the compile.
