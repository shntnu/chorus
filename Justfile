# chorus dev tasks. Run inside `nix develop` (or with direnv active).

# Default: list targets.
default:
    @just --list

# First-time setup: create conda env, install chorus, pull the Enformer oracle.
# Idempotent on re-runs (mamba/pip skip if satisfied; chorus setup resumes).
setup: env install oracle-enformer
    @echo "chorus ready. Try: just demo"

# Create the base chorus conda env from environment.yml.
env:
    micromamba create -y -n chorus -f environment.yml || \
        micromamba install -y -n chorus -f environment.yml

# Editable install of the chorus Python package into the env.
install:
    micromamba run -n chorus pip install -e .

# Pull the lightweight CPU-friendly starter oracle (~6 GB).
oracle-enformer:
    micromamba run -n chorus chorus setup --oracle enformer

# Add another oracle later. Usage: just oracle alphagenome
oracle name:
    micromamba run -n chorus chorus setup --oracle {{name}}

# Smoke test: predict wild-type DNase signal + scan all 3 alt alleles
# at chr11:5,247,500 (beta-globin locus, K562) using Enformer.
# Confirms the install actually works end-to-end.
demo:
    #!/usr/bin/env bash
    set -euo pipefail
    micromamba run -n chorus python - <<'PY'
    import chorus
    from chorus.utils import get_genome

    oracle = chorus.create_oracle(
        'enformer', use_environment=True,
        reference_fasta=str(get_genome('hg38')),
    )
    oracle.load_pretrained_model()

    wt = oracle.predict(('chr11', 5247000, 5248000), ['ENCFF413AHU'])
    print(f"WT DNase mean signal (K562, chr11:5,247,000-5,248,000): "
          f"{wt['ENCFF413AHU'].values.mean():.3f}")

    effects = oracle.predict_variant_effect(
        'chr11:5247000-5248000', 'chr11:5247500',
        ['C', 'A', 'G', 'T'], ['ENCFF413AHU'],
    )
    n_alts = len(effects['predictions']) - 1
    print(f"Scored {n_alts} alt alleles at chr11:5,247,500: "
          f"{list(effects['predictions'].keys())}")
    PY

# Drop into a python REPL inside the chorus env.
repl:
    micromamba run -n chorus python

# Nuke the conda env entirely. Models stay in HuggingFace cache.
clean:
    micromamba env remove -y -n chorus
