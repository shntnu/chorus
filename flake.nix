{
  description = "Chorus dev shell — micromamba-based; chorus owns the per-oracle conda envs";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    systems.url = "github:nix-systems/default";
    flake-utils.url = "github:numtide/flake-utils";
    flake-utils.inputs.systems.follows = "systems";
  };

  outputs =
    { self, nixpkgs, flake-utils, systems, ... }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          system = system;
          config.allowUnfree = true;
        };

        # System libs that pip/conda-installed binaries dlopen at runtime.
        # NIX_LD makes them findable to non-Nix-built executables.
        libList = [
          pkgs.zlib
          pkgs.stdenv.cc.cc
          pkgs.bzip2
          pkgs.xz
          pkgs.openssl
        ];
      in
      {
        devShells.default = pkgs.mkShell {
          packages = [
            pkgs.micromamba
            pkgs.just
            pkgs.git
            pkgs.git-lfs
            pkgs.which
            pkgs.curl
            pkgs.wget
          ] ++ libList;

          NIX_LD = pkgs.runCommand "ld.so" { } ''
            ln -s "$(cat '${pkgs.stdenv.cc}/nix-support/dynamic-linker')" $out
          '';
          NIX_LD_LIBRARY_PATH = pkgs.lib.makeLibraryPath libList;

          shellHook = ''
            # NVIDIA driver libs (libcuda.so) live here on NixOS — pip-installed
            # CUDA wheels (jax[cuda12], tensorflow) need this path. Pattern from
            # shntnu-neusis/templates/python-pixi/flake.nix:31.
            [ -d /run/opengl-driver/lib ] && \
              export LD_LIBRARY_PATH=/run/opengl-driver/lib:$LD_LIBRARY_PATH
            export LD_LIBRARY_PATH=$NIX_LD_LIBRARY_PATH:$LD_LIBRARY_PATH

            # Keep mamba state inside the repo so it doesn't pollute $HOME.
            export MAMBA_ROOT_PREFIX="$PWD/.mamba"
            eval "$(micromamba shell hook --shell bash)"

            if [ ! -d "$MAMBA_ROOT_PREFIX/envs/chorus" ]; then
              echo
              echo "chorus dev shell ready. First-time setup: run 'just setup'"
              echo "(see 'just' for other targets)"
              echo
            else
              micromamba activate chorus 2>/dev/null || true
              echo "chorus env active. Run 'just' to list targets, or 'chorus --help'."
            fi
          '';
        };
      }
    );
}
