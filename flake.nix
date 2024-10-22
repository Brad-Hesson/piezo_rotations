{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };
  outputs = flakes: flakes.flake-utils.lib.eachDefaultSystem (system:
    let
      pkgs = import flakes.nixpkgs {
        inherit system;
        config.allowUnfree = true;
      };
    in
    {
      devShell = with pkgs; mkShell {
        packages = [
          (pkgs.python3.withPackages (ps: [
            ps.matplotlib
            ps.numpy
            ps.scikit-image
            ps.plotly
          ]))
        ];
      };
    }
  );
}

