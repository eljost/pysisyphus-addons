final: prev: {
  python3 = prev.python3.override (old: {
    packageOverrides = prev.lib.composeExtensions (old.packageOverrides or (_: _: { })) (pfinal: pprev: {
      pysisyphus-addons = pfinal.callPackage ./pysisyphus-addons.nix { };
    });
  });
}
