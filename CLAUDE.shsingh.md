# CLAUDE.shsingh.md

Personal context for this fork. Not for upstream.
Read alongside upstream CLAUDE.md - that one is the project's; this
one is mine.

## Why this fork exists

This fork is exploratory work in service of the `pign-cdg` project at
<https://github.com/broadinstitute/pign-cdg> (private repo; clone as
sibling at `../pign-cdg` - see Operational notes). PIGN-CDG is a rare
congenital disorder of glycosylation; the project supports drug
repurposing and variant-impact research for the index patient and the
broader ~140-family cohort. Full project context including biology,
patient details, collaborators, and data inventory lives in
`docs/project-context.md` of that repo. Read it first if you need
background on PIGN-CDG itself.

The question driving this fork: can chorus contribute meaningfully to
that work, and if so where?

## Integration with pign-cdg - honest assessment

What chorus does well that could help:

- **PIGN regulatory landscape.** Predict PIGN gene expression across
  cell types and tissues (chorus has thousands of CAGE / RNA-seq
  tracks). Useful for grounding which cell types best model the
  affected tissue (neural / muscle for MCAHS1).
- **Non-coding variants in the extended cohort.** If any of the ~140
  PIGN-CDG families carry promoter / UTR / intronic VUS, chorus can
  predict regulatory effects across cell types. VarChAMP (imaging)
  and VISTA (yeast functional) don't cover this class.
- **Cross-modality demonstration.** Pick one variant where we have
  imaging (VarChAMP), functional (VISTA), and regulatory (chorus)
  signals. Demonstrate joint interpretation. This is the larger frame
  beyond pign-cdg specifically.

What chorus does NOT do that pign-cdg actually needs:

- The index patient's `p.Arg95Gln` is a coding missense. Chorus
  predicts regulatory effects, not protein-function effects. Use
  AlphaFold / ESM / VarChAMP imaging for that.
- The maternal `exon 15-25 deletion` is a structural variant. Chorus
  is built for SNVs.
- Drug-response prediction (Perlara yeast hits, ascorbyl palmitate
  etc.) is outside chorus's scope.

The most likely useful output of this exploration: a short writeup
saying which classes of pign-cdg questions chorus is and isn't useful
for, with a concrete worked example for the "useful" cases.

## Goals for the exploration

- Get chorus running end-to-end on real hardware. Smoke test with
  Enformer, then add AlphaGenome (the comprehensive one) to actually
  exercise the Linux GPU box.
- Run one realistic pign-cdg-shaped query: predict PIGN expression
  across the available cell-type tracks; compare against what's known
  about tissue-specific PIGN biology.
- If interesting, score the small set of regulatory-region PIGN
  variants from ClinVar that exist (most known PIGN pathogenic
  variants are coding, but check).
- Form an opinion on whether chorus belongs in the pign-cdg toolkit
  or is a parallel exploration.

## Operational notes

- **Sibling-clone pign-cdg.** Clone <https://github.com/broadinstitute/pign-cdg>
  at `../pign-cdg` relative to this repo, so paths like
  `../pign-cdg/docs/project-context.md` work for cross-referencing
  without absolute paths. Both repos then live side-by-side under
  the same parent directory.
- Targets: a Linux server with H100 GPUs (primary) and a Mac
  (smoke tests only).
- Env mgmt: chorus uses per-oracle conda envs because TF/PyTorch/JAX
  don't coexist. The Nix flake here provides micromamba + system libs;
  micromamba owns everything inside `.mamba/`.
- Setup: `direnv allow && just setup` for the Enformer-only fast path,
  or `chorus setup` (no flag) for all six oracles.
- HuggingFace token at `~/.cache/huggingface/token`; needed only for
  gated oracles (AlphaGenome).

## Logging convention

Day-by-day exploration notes in `LOG.shsingh.md` (oldest first,
matching the pign-cdg `docs/log.md` convention). Only promote to
`../pign-cdg/docs/log.md` when a finding is actionable for pign-cdg
itself - link back to the relevant `LOG.shsingh.md` entry or commit
hash so the trail stays connected.

## On the MCP question (parked)

Chorus ships an MCP server as its primary recommended interface
(Step 4 of the README). I'm not focused on this right now - the goal
is using chorus, not auditing its interface design. If a sharper
opinion forms after working with it, capture it here. Otherwise
skip.

## Posture toward upstream

- My commits to this fork's main: `flake.nix`, `flake.lock`,
  `Justfile`, `.envrc`, `.gitignore` additions, plus this file. None
  of these belong in an upstream PR by default - they're shaped to my
  workflow, not to upstream's.
- For any upstream contribution: cut a feature branch off
  `upstream/main` (not `origin/main`), cherry-pick only the relevant
  changes, open the PR from there.

## Open questions

- Is there a regulatory-region PIGN variant in the cohort worth
  a cherry-picked demonstration?
- What's the smallest joint result that would be worth showing
  the upstream maintainer?
