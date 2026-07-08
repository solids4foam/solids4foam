---
name: Release checklist
about: Track the steps for a new solids4foam release
title: "Release vX.Y checklist"
labels: release
assignees: ""
---

## Pre-release

- [ ] All target PRs merged into `development`; CI green (`buildAndTest`, `buildAndRegressionTest`)
- [ ] `./Alltest` and `./Alltest-regression` pass locally on each supported OpenFOAM/foam-extend version
- [ ] Update `CITATION.cff` (`version`, `date-released`, and DOI if a new one is minted)
- [ ] Update `CONTRIBUTORS.md` if new contributors since last release
- [ ] Tutorial/README docs reflect any new features or breaking changes

## Tag and publish

- [ ] Merge `development` into `master`
- [ ] Tag release on `master` (`vX.Y`) and push tag
- [ ] Create GitHub Release from tag with release notes (**not a draft** — publish directly to avoid the Zenodo draft-release breakage from a past release)
- [ ] [Optional] Create/update `CHANGELOG.md` in `master` branch with a summary of the main changes since the last release, where this file should ideally provide the full history, one section per version, newest first, back to v1.0
- [ ] Confirm Zenodo archive was created correctly from the GitHub Release (check the Zenodo record renders, metadata/DOI match `CITATION.cff`)

## Docker images

- [ ] Run `docker-release.yml` (workflow_dispatch) with the new version tag
- [ ] Pull each pushed image (`openfoam-v2412`, `openfoam-v2312`, `openfoam-9`, `foam-extend-4.1`) and sanity-check it runs a tutorial

## Website

- [ ] Confirm `update-submodule.yml` ran on push to `development`/`master` and solids4foam.github.io picked up the new commit
- [ ] Spot-check solids4foam.github.io renders correctly (no broken links/build)

## Announce

- [ ] Post release announcement (OpenFOAM/CFD forum, mailing list, etc., as applicable)
