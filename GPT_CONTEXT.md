# GPT Context

## Repository identity

This is the public `ANC-Research/IC-DMCANC-CPA` repository.

It contains a paper-specific MATLAB implementation of distributed multichannel ANC with intermittent communication and coprocessor-assisted data combination. Do not treat it as a general ANC library, an event-triggered asynchronous-communication project, feedback ANC, or output-constrained ANC.

Because this repository is public, do not add private repository names, unpublished methods, private data descriptions, internal results, credentials, or confidential future plans.

## Required reading before making changes

Before modifying code, read the latest version of:

1. `PROJECT_SCOPE.md`;
2. the root `README.md`;
3. `FedMCANC.m`;
4. the exact `FedDMCANC_case*.m` involved in the request;
5. directly affected comparison implementations;
6. the paper definition or existing implementation evidence relevant to an algorithm or architecture change.

Always use the latest requested branch. Do not rely only on previous conversations, cached code, or similarly named distributed ANC implementations from another repository.

## Implementation authority

- `FedMCANC.m` is the authoritative code for local WCFxLMS-style adaptation, compensation filters, MWD aggregation, ideal communication, common-period intermittent communication, and individual node communication periods.
- `FedDMCANC_case1.m` through `FedDMCANC_case5.m` are independent reproduction cases with different data and parameter purposes.
- `DMANC_CompensateSP.m` is the MGDFxLMS / gradient-transmission comparison used in this project.
- `McANC_FxLMS_SIMO.m` is the centralized 1 × 6 × 6 comparison used in this project.
- The committed `.mat` files define the current public simulation paths and compressor-noise asset.
- Keep unrelated cases and comparison methods read-only unless the request affects them.

## Task isolation

At the start of each task, identify:

- the target method: ideal communication, fixed built-in interval, common periodic communication, individual periodic communication, compensation-filter identification, local constraint, MWD fusion, coprocessor abstraction, baseline, data path, or plotting;
- the exact entry script and method name;
- the controlling sample rate, duration, node count, filter lengths, step sizes, `alpha`, communication interval, paths, and input data;
- whether the task changes equations, communication behaviour, hardware interpretation, dimensions, data dependencies, validation, or only presentation;
- which controller states, residual errors, communication events, plots, and saved variables are expected to change.

Do not assume that values or communication schedules from one case apply to another.

## Current structural assumptions

The current implementation is explicitly expanded for:

```text
1 reference × 6 secondary sources × 6 error sensors
```

Important shapes include:

```text
Wc:       (6, Lw + Lc - 1)
Wcsubopt: (6, Lw + Lc - 1)
Nabla:    (6, Lw + Lc - 1)
SecP:     (6, 6, Ls)
Dis:      (6, N)
yc:       (6, N)
C:        (6, 6, Lc)
```

Many controller, error, and fusion operations are expanded node by node. Do not claim arbitrary `K` support merely because the constructor accepts `node_num`. Generalization requires replacing and validating every six-node-specific operation and case assumption.

## Core algorithm invariants

Preserve these unless the user explicitly requests an algorithm change:

```text
e_k(n) = Dis_k(n) - y_k(n)
```

```text
W_k(n+1) = W_k(n)
           + mu * xf_k(n) * e_k(n)
           + alpha * mu * (Wcenter_k(n) - W_k(n))
```

```text
Nabla = Wc - Wcsubopt
```

Also preserve:

- sample-by-sample recursive local adaptation;
- full cross-secondary-path coupling in physical residual generation;
- diagonal secondary-path filtering for local updates;
- separate local-controller and center-point states;
- non-communication periods continuing local adaptation;
- compensation filtering of other-node weight differences before fusion;
- the selected tail segment after compensation filtering;
- local-controller reset after the relevant communication event;
- distinct ideal, common-period, and individual-period modes;
- equal experiment conditions when comparing communication settings.

Do not replace the local constraint, MWD fusion, or intermittent schedule with a centralized FxLMS update without explicitly defining a new method.

## Communication-mode rules

### Ideal mode

`FedMCANC_166_ideal` performs central aggregation inside every sample iteration. Treat it as the ideal-communication comparison, not an intermittent method.

### Built-in interval

`FedMCANC_166` currently uses:

```text
mod(i, 16000 * 10) == 0
```

The interval is tied to a hard-coded 16 kHz sample rate. A sample-rate or interval change must update the relationship explicitly.

### Common periodic communication

`FedMCANC_166_PCF(..., tc)` uses a shared period:

```text
mod(i, 16000 * tc) == 0
```

When changing this mode:

- state `tc` in seconds and samples;
- state the assumed sample rate;
- verify integer sample intervals;
- count communication rounds;
- compare convergence and final noise reduction with the ideal mode.

### Individual periodic communication

`FedMCANC_166_PCF_individual(..., tc)` assigns one period per node.

When changing this mode:

- verify `tc` has one valid value per node;
- identify which `Nabla(k,:)` values are refreshed at each event;
- explain how stale differences from non-communicating nodes are used;
- verify only the intended node's center point and local controller are reset;
- report per-node and total communication counts;
- compare heterogeneous and common-period results under comparable traffic.

Do not call a schedule asynchronous or heterogeneous without stating exactly whether communication times differ across nodes and how stale information is handled.

## Local constraint and fusion changes

When changing `alpha`, `Wcsubopt`, `Nabla`, compensation filters, or aggregation:

1. map every mathematical term to the exact MATLAB property and update line;
2. distinguish local adaptation from central aggregation;
3. verify lengths `Lw`, `Lc`, and `Lw + Lc - 1`;
4. verify sign and orientation of compensation filters;
5. check the tail segment selected after filtering;
6. state which nodes transmit and which nodes are updated;
7. preserve or explicitly redefine the center-point reset step;
8. compare controller norms, finite values, communication traffic, and residual error;
9. update all communication modes only when mathematically intended.

A center-point constraint or WCFxLMS term must not be described as a universal formal stability proof for every network, plant, and communication schedule.

## Coprocessor interpretation

The code models coprocessor assistance as a centralized MATLAB aggregation operation.

Do not claim that this repository contains:

- coprocessor firmware;
- real-time task scheduling;
- network drivers or protocols;
- hardware communication interfaces;
- fixed-point or embedded deployment;
- measured execution-time or processor-load results.

When modifying architecture documentation, distinguish clearly between:

- algorithmic central aggregation implemented in MATLAB;
- the real-time coprocessor architecture described in the paper;
- hardware code that is not present in this repository.

## Experiment-case boundaries

- `FedDMCANC_case1.m`: synthetic noise, ideal-network comparison.
- `FedDMCANC_case2.m`: included compressor noise, ideal-network comparison.
- `FedDMCANC_case3.m`: effect of common communication intervals.
- `FedDMCANC_case4.m`: effect of `alpha` at a fixed interval.
- `FedDMCANC_case5.m`: node-specific communication intervals.

Do not reuse results, parameters, or plots across cases without checking their input data, duration, step size, communication mode, and intended paper figure.

## Data and path rules

- The committed path assets and all five case scripts use the repository-root-relative directory `simulation path/`.
- Run the cases from the repository root, or ensure the repository root is available on the MATLAB path.
- Do not silently duplicate, rename, or move `.mat` assets.
- A future directory change must update all affected scripts and documentation together and then be validated from the repository root.
- Preserve the included compressor recording and simulation paths unless the user explicitly requests an asset change.
- Do not commit private recordings, measured paths, credentials, or machine-specific absolute paths.
- Treat generated plots and saved workspaces as experiment artifacts, not source-code authority.
- Record random seeds when reproducibility is introduced or changed.

## MATLAB change discipline

- Make the smallest coherent change that satisfies the request.
- Avoid broad stylistic rewrites of the expanded six-node reference code during an algorithm fix.
- Preserve class and method interfaces unless an interface change is required.
- Check row/column orientation and `reshape` order explicitly.
- Check `mod` conditions for integer sample intervals and sample-rate dependence.
- Keep algorithm, communication, path, and plotting changes separate where practical.
- Report every modified, added, renamed, or deleted file.
- Add a concise entry to `CHANGELOG.md` for significant user-visible changes.
- Update the root README when setup, paths, communication modes, dependencies, architecture interpretation, or reproduction instructions change.

## Cross-repository restrictions

Do not import code, parameters, or conclusions from other ANC repositories unless the user explicitly requests a documented dependency or comparison.

In particular, do not merge this implementation with the event-triggered asynchronous-communication project. The repositories represent different communication architectures and must remain independently traceable.

Similarly named WCFxLMS, MWD, FxLMS, compensation-filter, or communication routines in different repositories are independent implementations until their equations, schedules, and conventions have been compared.

## Verification

Use the strongest verification available for the change.

For static inspection, check:

- MATLAB syntax, class names, method names, and call signatures;
- actual asset paths and paths used by scripts;
- controller, center-point, difference, compensation-filter, and history dimensions;
- initialization of `Wc`, `Wcsubopt`, `Nabla`, and compensation filters;
- sample-rate dependence and integer validity of communication intervals;
- finite values in controllers, errors, compensation filters, and fused differences;
- continued local adaptation during non-communication intervals;
- node-specific updates affecting only intended nodes;
- comparison methods receiving the same disturbance, reference, paths, and step sizes.

For MATLAB execution, run the directly affected case first, for example:

```matlab
FedDMCANC_case1
FedDMCANC_case3
FedDMCANC_case5
```

As applicable, inspect:

- residual-noise curves for each node and their mean;
- ideal, intermittent, centralized, and MGDFxLMS comparisons;
- communication rounds per node and total communication traffic;
- local and center-point controller norms;
- compensation-filter identification errors;
- convergence and final noise-reduction trade-offs;
- results for different `alpha` and `Tc` values;
- regenerated plots and saved variables.

MATLAB toolbox availability and long sample-by-sample simulations may prevent execution in the current environment. Clearly separate code inspection from MATLAB results that must be produced locally.

## Response requirements

For every completed modification, provide:

- files changed;
- the purpose of each change;
- the affected communication mode and experiment case;
- original and updated equation or schedule mapping, when applicable;
- parameter, path, and dimension impact;
- coprocessor/hardware interpretation impact, when applicable;
- verification performed and results;
- MATLAB or local-runtime checks still required.