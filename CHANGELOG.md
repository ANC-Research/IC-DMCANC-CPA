# Changelog

This file records significant repository-level changes from its introduction onward. Earlier development history remains available in Git history and is not reconstructed here.

## Unreleased

### Added

- Added `PROJECT_SCOPE.md` to define the repository as a paper-specific intermittent-communication Distributed ANC project, identify the authoritative MATLAB implementation and experiment cases, document the six-node system, WCFxLMS-style local update, MWD fusion, communication modes, coprocessor abstraction, boundaries, and validation expectations.
- Added `GPT_CONTEXT.md` to define required reading, task isolation, six-node structural assumptions, communication-schedule rules, local-constraint and fusion boundaries, coprocessor interpretation, experiment-case separation, public-content restrictions, path handling, MATLAB change discipline, and verification requirements.
- Added this changelog for future user-visible algorithm, communication schedule, fusion, parameter, data-path, architecture, experiment, interface, and documentation changes.

### Fixed

- Corrected the acoustic-path load statements in `FedDMCANC_case1.m` through `FedDMCANC_case5.m` to use the committed repository directory `simulation path/` instead of the nonexistent `path/simulation path/` prefix.
- Updated `PROJECT_SCOPE.md` and `GPT_CONTEXT.md` to reflect the corrected repository-root-relative paths.