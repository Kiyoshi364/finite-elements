{
  description = "A flake with builds and shells for computacao-algebrica";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-24.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
        defaultShellHook = { shellName, ... }: ''
            export PS1="\n\[\033[1;32m\][${shellName}: \w]\n\$\[\033[0m\] "
            export PURE="$([[ $IN_NIX_SHELL = 'pure' ]] && echo 1 || echo 0)"
            echo "PURE=$PURE"
            echo -n '>> Welcome to ${shellName}! <<'
          '';
        myShell = input@{
            name,
            buildInputs ? [],
            toShellName ? n: "${n}-shell",
            shellHook ? defaultShellHook,
          }: pkgs.mkShell {
            name = toShellName name;

            buildInputs = buildInputs;

            shellHook = shellHook
              (input // {
                inherit defaultShellHook;
                shellName = toShellName name;
              });
          };
      in {
        devShells = {
          default = self.devShells.${system}.julia;

          julia = myShell {
            name = "julia";
            buildInputs = with pkgs; [
              julia
            ];
          };

          pdf = myShell {
            name = "zathura";
            buildInputs = with pkgs; [
              zathura
            ];
          };

          tex = myShell {
            name = "tex";
            buildInputs = with pkgs; [
              (texlive.combine
                { inherit (texlive)
                  scheme-small
                ; }
              )

              (aspellWithDicts (dicts: with dicts; [
                pt_BR
              ]))

              poppler_utils
              zathura
            ];
          };

          tex-summary = myShell {
            name = "tex-summary";
            buildInputs = with pkgs; [
              (texlive.combine
                { inherit (texlive)
                  scheme-small

                  biblatex
                  csquotes

                  catchfile
                  transparent
                  trimspaces
                  svg
                ; }
              )

              (aspellWithDicts (dicts: with dicts; [
                en
                pt_BR
              ]))

              inkscape

              poppler_utils
              zathura
            ];
          };

          typst = myShell {
            name = "typst";
            buildInputs = with pkgs; [
              typst

              (aspellWithDicts (dicts: with dicts; [
                pt_BR
              ]))

              poppler_utils
              zathura
            ];
          };

          md = myShell {
            name = "markdown";
            buildInputs = with pkgs; [
              marker
            ];
          };
        };
      });
}
