const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaCorr(alphaControls.get<label>("nAlphaCorr"));

label nAlphaSubCycles(alphaControls.get<label>("nAlphaSubCycles"));

bool MULESCorr(alphaControls.getOrDefault("MULESCorr", false));

// Apply the compression correction from the previous iteration
// Improves efficiency for steady-simulations but can only be applied
// once the alpha field is reasonably steady, i.e. fully developed
bool alphaApplyPrevCorr
(
    alphaControls.getOrDefault("alphaApplyPrevCorr", false)
);

// Isotropic compression coefficient
scalar icAlpha
(
    alphaControls.getOrDefault<scalar>("icAlpha", 0)
);

// Shear compression coefficient
scalar scAlpha
(
    alphaControls.getOrDefault<scalar>("scAlpha", 0)
);
