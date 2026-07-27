# Project Scope

## Purpose

`IC-DMCANC-CPA` is the public MATLAB reference repository for **Distributed Multichannel Active Noise Control with Intermittent Communication and Coprocessor-Assisted Data Combination**.

The repository accompanies the paper:

> Implementation of distributed multichannel active noise control with intermittent communication and coprocessor assisted data combination

This is a paper-specific Distributed ANC project. It is not a general ANC algorithm collection and must remain independently traceable from other communication-efficient ANC repositories.

## Authoritative implementation

The proposed controller and communication methods are implemented in:

```text
FedMCANC.m
```

The reproduction cases are:

```text
FedDMCANC_case1.m
FedDMCANC_case2.m
FedDMCANC_case3.m
FedDMCANC_case4.m
FedDMCANC_case5.m
```

Comparison implementations are:

```text
DMANC_CompensateSP.m
McANC_FxLMS_SIMO.m
```

The MATLAB source and the associated paper are the authority for algorithm definitions, communication schedules, experiment conditions, and result interpretation.

## Current system model

The included implementation is explicitly expanded for:

```text
1 reference × 6 secondary sources × 6 error sensors
```

Important internal shapes are:

```text
Wc:       6 × (Lw + Lc - 1)
Wcsubopt: 6 × (Lw + Lc - 1)
Nabla:    6 × (Lw + Lc - 1)
SecP:     6 × 6 × Ls
Dis:      6 × N
yc:       6 × N
C:        6 × 6 × Lc
```

The current experiment scripts generally use:

```text
Fs       = 16000 Hz
wLen     = 512
sLen     = 256
Numnode  = 6
cLen     = 33
muc      = 1e-5
```

The controller step size, duration, constraint parameter `alpha`, noise source, and communication interval vary by case and must be read from the selected script.

Although the constructor accepts `node_num`, many controller and fusion operations are manually expanded for six nodes. Arbitrary-node support must not be claimed without generalizing and validating every hard-coded operation.

## Repository paths

The committed acoustic-path assets are located at:

```text
simulation path/PrimaryPath_1x6.mat
simulation path/SecondaryPath_6x6.mat
```

All five `FedDMCANC_case*.m` scripts use these repository-root-relative paths. The cases are intended to be launched with the repository root as the MATLAB working directory, or with the repository root available on the MATLAB path.

Do not move or duplicate the `.mat` assets silently. A future directory change must update all affected scripts and documentation together.

## Algorithm components

### Compensation-filter identification

`FedMCANC.CompensateSP` identifies compensation filters for cross-secondary-path contributions.

For each off-diagonal secondary path, white noise is passed through the cross path. A Filtered-x LMS system using the corresponding diagonal path produces a compensation filter stored in:

```text
C(m, k, :)
```

The implementation applies a sign reversal and coefficient flip before storing the result. Any change to the identification signal, diagonal/cross-path mapping, filter orientation, sign, length, or convergence procedure changes the fusion implementation.

### Local WCFxLMS-style adaptation

Each node performs local sample-by-sample adaptation using its local error and diagonal secondary path. The implemented update has the form:

```text
W_k(n+1) = W_k(n)
           + mu * xf_k(n) * e_k(n)
           + alpha * mu * (Wcenter_k(n) - W_k(n))
```

where:

- `Wc` is the current local controller state;
- `Wcsubopt` is the center-point or constrained reference state;
- `alpha` controls attraction to the center point;
- `xf_k` is formed using the diagonal secondary path;
- the physical residual error uses the complete 6 × 6 secondary-path plant.

Changing the sign, scaling, filtered-reference definition, center-point role, or interpretation of `alpha` is an algorithm change.

### Weight-difference and MWD fusion

At a communication event, the implementation forms:

```text
Nabla = Wc - Wcsubopt
```

The central aggregation stage combines each node's own weight difference with compensation-filtered differences from the other nodes. The center-point controller is updated, then selected local controllers are reset to the new center point.

This is the code-level representation of Mixed Weight Difference data combination. Changes to the transmitted quantity, compensation filtering, selected tail segment, aggregation order, node participation, or reset behaviour are communication-algorithm changes.

### Coprocessor-assisted architecture

The MATLAB code represents the coprocessor role through a centralized aggregation block that performs weight-difference fusion.

The repository does **not** currently contain coprocessor firmware, hardware interfaces, scheduling code, fixed-point implementation, communication drivers, or real-time deployment software. Hardware and real-time architecture details are described in the paper and must not be presented as executable code contained here.

## Communication modes

### Ideal-network comparison

```text
FedMCANC_166_ideal
```

This method performs the central aggregation during every sample iteration and serves as the ideal-communication comparison.

### Fixed built-in intermittent communication

```text
FedMCANC_166
```

The current method performs central aggregation when:

```text
mod(i, 16000 * 10) == 0
```

which corresponds to a 10-second interval only when the sample rate is 16 kHz. This hard-coded sample-rate dependence must be documented if changed.

### Parameterized common communication frequency

```text
FedMCANC_166_PCF(..., tc)
```

All nodes communicate at the same period using:

```text
mod(i, 16000 * tc) == 0
```

The value `tc` is expressed in seconds only under the current 16 kHz assumption.

### Individual node communication frequencies

```text
FedMCANC_166_PCF_individual(..., tc)
```

Each node has its own communication period `tc(k)`. At a node-specific event, the implementation updates that node's transmitted difference, center point, and local controller while retaining the most recently available differences from other nodes.

This method represents heterogeneous intermittent communication. Changes to stale-difference handling, synchronization, event ordering, or which node state is reset must be documented explicitly.

## Experiment cases

| Case | Purpose |
|---|---|
| `FedDMCANC_case1.m` | Ideal-network comparison using synthetic band-limited random noise and simulation paths. |
| `FedDMCANC_case2.m` | Ideal-network comparison using the included compressor recording. |
| `FedDMCANC_case3.m` | Effect of a common communication interval, including `Tc = [0.1, 0.5, 1, 3]` seconds in the current script. |
| `FedDMCANC_case4.m` | Effect of the WCFxLMS center-point constraint parameter `alpha` at a fixed communication interval. |
| `FedDMCANC_case5.m` | Heterogeneous node-specific communication intervals, currently `Tc = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]` seconds. |

Case parameters and saved outputs are independent experiment definitions. Do not copy values between cases without checking the intended comparison.

## Comparison methods

| File | Role |
|---|---|
| `FedMCANC.m` | Proposed intermittent-communication methods, ideal-network comparison, local WCFxLMS-style adaptation, and central MWD fusion. |
| `DMANC_CompensateSP.m` | MGDFxLMS / gradient-transmission distributed comparison with compensation filters. |
| `McANC_FxLMS_SIMO.m` | Centralized 1 × 6 × 6 multichannel FxLMS comparison. |
| `FedDMCANC_case*.m` | Paper reproduction scripts and parameter definitions. |

The comparison implementations are authoritative only for this repository's reproduction and must not be treated as general-purpose cross-repository libraries.

## In scope

- Intermittent-communication distributed multichannel ANC.
- WCFxLMS-style local constrained adaptation.
- Common and node-specific communication periods.
- Mixed Weight Difference exchange and compensation-filtered fusion.
- Software modelling of coprocessor-assisted central aggregation.
- Ideal-network, centralized, and MGDFxLMS comparisons.
- Six-node simulation paths and the included compressor-noise experiment.
- Communication-frequency, constraint-parameter, residual-error, and noise-reduction evaluation.
- MATLAB documentation and clearly scoped bug fixes.

## Out of scope

- General-purpose ANC algorithm collections.
- Event-triggered asynchronous communication methods maintained as a separate project.
- Feedback ANC or output-constrained MOV-FxLMS research.
- Claims of arbitrary node-count support without code generalization and validation.
- Coprocessor firmware, communication drivers, packet protocols, fixed-point deployment, or real-time hardware code not present in the repository.
- Packet loss, transmission delay, quantization, topology changes, or network scheduling unless added as a documented extension.
- Python, C/C++, embedded, or firmware ports unless introduced as a separate implementation track.
- Private acoustic data, credentials, machine-specific absolute paths, or confidential results.

## Core invariants

Unless the user explicitly requests an algorithm change, preserve:

- residual-error convention `e = Dis - y`;
- sample-by-sample recursive local adaptation;
- full 6 × 6 physical secondary-path coupling in residual generation;
- diagonal secondary-path filtering for each local update;
- separate local-controller and center-point states;
- WCFxLMS center attraction controlled by `alpha`;
- weight-difference definition `Nabla = Wc - Wcsubopt`;
- compensation-filtered fusion of other-node differences;
- reset of the intended local controller after communication;
- distinct ideal, common-period, and individual-period communication modes;
- equal paths, inputs, step sizes, and evaluation windows when comparing methods.

A center-point constraint or intermittent communication schedule must not be described as a general formal proof of closed-loop stability.

## MATLAB dependencies

The current code uses MATLAB functions that may require relevant toolboxes, including:

- `dsp.FilteredXLMSFilter`;
- `awgn`;
- `fir1`, `filter`, and `smooth`;
- plotting and MAT-file loading functions.

Toolbox-dependent execution results must be separated from static code inspection when the required MATLAB environment is unavailable.

## Validation expectations

Changes should use the strongest applicable checks:

1. inspect `FedMCANC.m` and the affected `FedDMCANC_case*.m` together;
2. identify whether the change affects local adaptation, compensation filters, ideal communication, common-period communication, node-specific communication, central aggregation, or a baseline;
3. verify dimensions of `Wc`, `Wcsubopt`, `Nabla`, `SecP`, `C`, `Dis`, `yc`, and `e`;
4. verify the actual committed data path and the paths used by scripts;
5. check sample-rate dependence of all communication intervals;
6. verify finite controller weights, errors, compensation filters, and fused differences;
7. verify that non-communication intervals continue local adaptation;
8. verify that individual-period communication updates only the intended node at each event;
9. compare communication counts, convergence, final residual noise, and performance against the ideal-network case;
10. rerun the directly affected MATLAB case and regenerate plots when MATLAB is available;
11. document changes to equations, communication schedules, `alpha`, data paths, experiment assets, or result generation;
12. report every modified, added, renamed, or deleted file.

When MATLAB execution is unavailable, clearly separate code-level validation from reproduction results that still require local MATLAB and toolbox execution.