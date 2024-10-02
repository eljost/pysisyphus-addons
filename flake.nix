{
  description = "Compiled addons for pysisyphus";

  inputs = {
    flake-utils.url = "github:numtide/flake-utils";
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    sympleints.url = "github:eljost/sympleints";
  };

  outputs = { self, nixpkgs, sympleints, flake-utils }:
    flake-utils.lib.eachDefaultSystem
      (system:
        let
          overlay = import ./nix/overlay.nix;
          pkgs = import nixpkgs {
            inherit system;
            overlays = [ sympleints.overlays.default overlay ];
          };
        in
        {
          packages.default = pkgs.python3.pkgs.pysisyphus-addons;

          formatter = pkgs.nixpkgs-fmt;
        }) // {
      overlays.default = import ./nix/overlay.nix;
    };
}
