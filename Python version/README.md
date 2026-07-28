# IC-DMCANC-CPA Python reference

This directory contains an independent Python reference implementation of the
six-node simulations in the repository. The existing MATLAB files remain the
authoritative definitions. This port is intended for reproducible algorithm
inspection and offline simulation; it is not a real-time, embedded, firmware,
or coprocessor implementation.

## Scope

The first release is deliberately fixed at:

```text
1 reference × 6 secondary sources × 6 error sensors
NUM_NODES = 6
DTYPE = numpy.float64
```

The implementation does **not** claim support for arbitrary node counts. It
does not add GPU, PyTorch, JAX, block/frequency-domain FxLMS, real-time audio,
fixed point, packet loss, delay, quantization, or event-triggered
communication. Those are possible future tracks, not features of this strict
reference.

## MATLAB-to-Python map

| MATLAB authority | Python implementation |
|---|---|
| `FedMCANC.CompensateSP` | `ic_dmcanc/compensation_filter.py` |
| `FedMCANC_166_ideal` | `run_fed_mcanc(..., mode="ideal")` |
| `FedMCANC_166` | `mode="fixed_10_seconds"` |
| `FedMCANC_166_PCF` | `mode="common_periodic"` |
| `FedMCANC_166_PCF_individual` | `mode="individual_periodic"` |
| `McFxLMS_SIMO_166` | `ic_dmcanc/centralized_fxlms.py` |
| `DMANC_gradient_166` | `ic_dmcanc/mgdfxlms.py` |
| `FedDMCANC_case1.m` … `case5.m` | `cases/case1.py` … `case5.py` |

Repository assets are loaded in place; they are not moved or duplicated:

```text
../simulation path/PrimaryPath_1x6.mat
../simulation path/SecondaryPath_6x6.mat
../compressor_16kHz.mat
```

The loaders require the actual committed variables and shapes:

```text
Primary_path   (6, 512)
Secondary_path (6, 6, 256)
tr             (279378, 1)
yr             (279378, 1)
```

No unconditional `squeeze` or machine-specific absolute path is used.

## Core algorithm

For source \(k\), error sensor \(m\), and the complete physical secondary-path
plant:

\[
y'_m(n)=\sum_{k=1}^{6} S_{m,k} * y_k(n)
\]

\[
e_m(n)=d_m(n)-y'_m(n)
\]

Each local controller uses only its corresponding diagonal path to form the
filtered reference:

\[
W_k(n+1)=W_k(n)+\mu_w x_{f,k}(n)e_k(n)
+\alpha\mu_w\left(W_{\mathrm{center},k}(n)-W_k(n)\right).
\]

The Python sample order matches `FedMCANC.m`:

1. insert the current reference sample into the controller buffer;
2. generate all six control outputs from the current local weights;
3. pass those outputs through the complete 6×6 secondary plant;
4. compute `residual_error = disturbance - secondary_output`;
5. filter the current reference through each diagonal secondary path;
6. update all six local controllers with the positive-sign rule above;
7. test communication using the MATLAB sample index `n + 1`;
8. aggregate/reset the applicable controller state;
9. record end-of-sample diagnostics.

`local_weights` and `center_weights` are independent arrays with shape:

```text
(6, w_len + c_len - 1)
```

Only the first `w_len` coefficients generate the acoustic control output, as in
MATLAB. The longer state is required for filtered-reference updates and MWD
fusion.

## Compensation-filter identification

Each call reproduces the structure of the two MATLAB `CompensateSP` methods:

- one white-noise sequence of length `200_000` by default;
- the same sequence is shared by all 30 off-diagonal path pairs;
- every pair starts with independent zero adaptive-filter state;
- `m == k` is skipped;
- the cross path produces the system-identification desired signal;
- the corresponding diagonal path is both the physical and estimated
  secondary path for Filtered-x LMS;
- the sample error is `desired - secondary_output`;
- coefficients update sample by sample with the filtered reference.

The raw identified coefficients are transformed exactly once:

```python
compensation_filter = -identified_coefficients[::-1]
```

The negative sign aligns system-identification output with the ANC/MWD sign
convention. The reversal matches the later filtering and final-tail
selection. Fusion must not reverse the coefficients again.

Cases 1 and 2 call compensation identification separately for MGDFxLMS and
FedDMCANC, as the MATLAB scripts do. Python uses deterministic independent
streams:

```text
MGDFxLMS: comp_id_seed
FedDMCANC: comp_id_seed + 1
```

Within each call, all off-diagonal pairs still share one white-noise vector.
Cases 3–5 require only the FedDMCANC identification and use `comp_id_seed`.

## MWD fusion and filter direction

At a common communication event:

```python
weight_difference = local_weights - center_weights
```

For an other-node contribution, the literal reference operation is:

```python
a = scipy.signal.lfilter(c, [1.0], weight_difference)
a_tail = a[-w_len:]
```

The required dimensions are:

```text
len(weight_difference) = w_len + c_len - 1
len(c)                 = c_len
len(a_tail)            = w_len
```

The accelerated loop evaluates those exact final `w_len` FIR outputs directly.
`tests/test_compensation_filter.py` checks it against `lfilter` tap by tap. It
does not use `convolve(..., mode="same")` or center cropping.

Each target center receives its own first-`w_len` difference, followed by the
compensation-filtered tail from every other node. Common/ideal aggregation
uses a full six-node difference snapshot before any center update.

## Communication modes

| Python mode | MATLAB method | Trigger |
|---|---|---|
| `ideal` | `FedMCANC_166_ideal` | every sample |
| `fixed_10_seconds` | `FedMCANC_166` | `fs * 10` samples |
| `common_periodic` | `FedMCANC_166_PCF` | `fs * tc_seconds` |
| `individual_periodic` | `FedMCANC_166_PCF_individual` | one period per node |

MATLAB uses one-based indexing:

```matlab
for i = 1:N
    if mod(i, interval_samples) == 0
```

Python therefore uses:

```python
matlab_sample_index = python_sample_index + 1
```

There is no event at Python sample zero unless the MATLAB sample index itself
reaches the interval.

`COMM_INTERVAL_ROUNDING = "require_integer"` is enforced. Decimal seconds are
converted with exact rational arithmetic. A non-integer `fs * tc_seconds`
raises an error; it is never silently floored, ceiled, rounded, or truncated.

In individual mode:

- only the triggering node refreshes its transmitted difference;
- all non-triggering differences remain stale;
- triggers at one sample execute in node order `0, 1, ..., 5`;
- each trigger immediately updates that node's center and resets only that
  local controller;
- a later trigger at the same sample sees differences refreshed by earlier
  nodes.

The event CSV uses one row per ideal/common aggregation round (`node=all`) and
one row per node trigger in individual mode. `summary.json` records both
per-node upload counts and the total number of node uploads.

## Cases

| Case | Input and purpose | MATLAB defaults |
|---|---|---|
| 1 | synthetic 100–1000 Hz; ideal centralized/MGDFxLMS/FedDMCANC comparison | 90 s, `mu_w=1e-6`, `alpha=1000` |
| 2 | first 14 s of `yr`; ideal baseline comparison | 14 s, `mu_w=3e-6`, `alpha=1000` |
| 3 | common communication-period sweep | `Tc=(0.1, 0.5, 1, 3)` s |
| 4 | center-attraction sweep at `Tc=0.5` s | `alpha=(300,600,1000,2000,5000,2.1e6)` |
| 5 | heterogeneous periods | `Tc=(0.2,0.3,0.4,0.5,0.6,0.7)` s |

Synthetic cases retain MATLAB's `0:1/Fs:T` length of `T*Fs + 1` samples. Their
display/evaluation window remains the first `T*Fs` samples. Case 2 uses exactly
`yr(1:T*Fs)`.

Case 3 and Case 4 use:

```text
SWEEP_STATE_MODE = "independent"
```

Every `Tc` or `alpha` starts from the same zero local/center state and the same
identified compensation-filter set. A chained state from a preceding sweep
point is not supported. The current committed MATLAB scripts call the base
`FedDMCANC` object for each point; this Python policy makes that intended
independence explicit.

## Comparison algorithms

Implemented:

- the current 1×6×6 path of `McFxLMS_SIMO_166`;
- the current 1×6×6 path of `DMANC_gradient_166`;
- the proposed four FedDMCANC communication modes.

Not implemented:

- historical 1×2×2 and 1×4×4 expansions not used by the five cases;
- fixed-controller/initial-controller variants not used by the cases;
- variable-delay MGDFxLMS variants not used by the cases;
- DFxLMS and ADFxLMS, which are external baselines in the MATLAB comments.

DFxLMS/ADFxLMS are not executed, plotted, or replaced with placeholder data.
An external result can be integrated later without changing the core
algorithms:

```python
# External baseline interface, implementation not included:
# diffusion_error, augmented_diffusion_error = ...
```

## Installation

Python 3.10 or newer is required.

```bash
cd "Python version"
python -m pip install -r requirements.txt
```

For full-duration simulations, install the optional Numba acceleration:

```bash
python -m pip install ".[fast]"
```

Tests can run with the standard library's `unittest`. Optional pytest support:

```bash
python -m pip install ".[test]"
```

## Running

From `Python version/`:

```bash
python -m cases.case1
python -m cases.case2
python -m cases.case3
python -m cases.case4
python -m cases.case5
```

Direct scripts are also supported:

```bash
python cases/case1.py
```

From the repository root:

```bash
python "Python version/cases/case1.py"
```

A short end-to-end check uses the real MAT assets but reduced adaptive lengths
and sample counts:

```bash
python -m cases.case1 --quick
```

Inspect all options:

```bash
python -m cases.case3 --help
```

Example overrides:

```bash
python -m cases.case3 --tc 0.1 0.5 1 3 --no-save-control
python -m cases.case4 --tc 0.5 --alpha-values 300 600 1000 2000
python -m cases.case5 --tc 0.2 0.3 0.4 0.5 0.6 0.7
```

## Main parameters

| Configuration field | Default | Meaning and valid values |
|---|---:|---|
| `fs` | 16000 | Sampling rate. The committed paths/recording are 16 kHz. |
| `duration_seconds` | 90 | Synthetic duration; case 2 overrides to 14 s. `fs*T` must be integer. |
| `num_nodes` | 6 | Must remain exactly 6. |
| `w_len` | 512 | Acoustic controller FIR taps. |
| `s_len` | 256 | Secondary path taps; must match the MAT asset. |
| `c_len` | 33 | Compensation-filter taps. |
| `mu_w` | `1e-6` | Nonnegative controller step size; case 2 uses `3e-6`. |
| `mu_c` | `1e-5` | Nonnegative compensation-ID step size. |
| `alpha` | 1000 | Nonnegative local center-attraction coefficient. |
| `communication_mode` | `ideal` | One of the four modes above. |
| `tc_seconds` | case-specific | One common value or exactly six individual values. |
| `random_seed` | 0 | Synthetic white-noise seed. |
| `awgn_seed` | 1 | Measured-SNR AWGN seed. |
| `comp_id_seed` | 0 | Compensation-ID base seed. |
| `awgn_snr_db` | 40 | AWGN uses the actual mean-square signal power without de-meaning. |
| `synthetic_low_hz` | 100 | Synthetic band lower edge. |
| `synthetic_high_hz` | 1000 | Synthetic band upper edge. |
| `synthetic_fir_order` | 63 | `fir1(63,...)` maps to 64 SciPy `firwin` taps. |
| `comp_id_num_samples` | 200000 | White-noise samples in each `CompensateSP` call. |
| `comp_id_convergence_stride` | 1000 | Saved system-ID diagnostic spacing only. |
| `diagnostic_stride` | 1600 | Saved controller/center norm spacing only. |
| `use_numba` | true | Use Numba when installed; formulas and state order are unchanged. |

Output switches are in `OutputConfig` and have matching CLI flags:

```text
save_full_error              True
save_full_control_output     True
save_full_secondary_output   False
save_full_comp_id_error      False
save_communication_events    True
save_controller_history      False
save_plots                   True
compress_npz                 True
```

## Outputs

Each entry point writes `outputs/caseN/` by default:

```text
config.json
summary.json
results.npz
metrics.npz
communication_events__*.csv
residual_level_per_node.png
residual_level_mean.png
controller_norms.png
compensation_identification.png
```

`results.npz` contains the clean/noisy reference, generated AWGN, disturbance,
paths, compensation-ID white noise, raw coefficients, final `C`, full residual
and optional control arrays, final weights, differences, communication arrays,
and norm diagnostics. Cases 1 and 2 save the independent MGDFxLMS and
FedDMCANC compensation results under distinct keys.

`metrics.npz` retains both unsmoothed per-sample relative levels and
MATLAB-style moving-average display metrics. `smooth`/`matlab_smooth` is never
used inside a controller update.

`summary.json` includes:

- finite-value checks;
- per-node and average final-window MSE/noise reduction;
- final controller and center norms;
- engine and runtime;
- communication periods, per-node counts, total uploads, and first/last event;
- input shapes and random seeds;
- estimated uncompressed result-array memory.

## Tests and MATLAB fixture

Run all local tests:

```bash
python -m unittest discover -s tests -v
```

The tests cover shape checks, `fir1`/AWGN/smoothing behavior, compensation-ID
update order, `-flip`, literal MWD tail selection, communication indexing and
counts, stale individual differences, same-sample node order, local update
sign, full plant coupling, independent initialization, both comparison
algorithms, and an independent MATLAB-style sample loop.

For a true MATLAB/Python short-trace comparison, first run in MATLAB:

```matlab
cd("Python version/matlab")
export_equivalence_fixture("matlab_equivalence.mat")
```

Then from `Python version/` on Linux/macOS:

```bash
IC_DMCANC_MATLAB_FIXTURE=matlab/matlab_equivalence.mat \
python -m unittest tests.test_matlab_equivalence -v
```

PowerShell:

```powershell
$env:IC_DMCANC_MATLAB_FIXTURE = "matlab\matlab_equivalence.mat"
python -m unittest tests.test_matlab_equivalence -v
```

The fixture contains fixed reference, disturbance, paths, initial state,
control output, residual error, and final controller states.

## Performance and memory

The readable kernels use the exact recursive order and are suitable for short
verification. A 90-second, 16 kHz, 512-tap simulation contains 1,440,001
samples and should use the optional Numba path. Block updates are intentionally
not used.

One `6 × 1,440,001` float64 signal is about 69 MiB. Saving disturbance,
residual, and control already requires about 207 MiB before paths, metrics, and
additional algorithms. Case 1 with three algorithms can require several
hundred MiB. Use `--no-save-control`, `--no-save-error`, or separate case runs
when memory is limited. Full `6 × 6 × 200000` system-ID error is disabled by
default.

## Validation status and known differences

Completed in the development environment:

- static variable/shape and update-order checks;
- 25 unit tests: 23 passed, 2 optional tests skipped;
- ideal, common-periodic, and individual-periodic short traces matched an
  independent MATLAB-style Python loop;
- centralized, MGDFxLMS, and FedDMCANC short runs were finite;
- all five `--quick` case entry points completed with finite algorithm
  outputs; the Case 1 end-to-end run also produced JSON, NPZ, CSV, and plots.

Not completed in the development environment:

- MATLAB execution was unavailable, so the exported MATLAB fixture test was
  skipped;
- Numba was unavailable, so compiled/reference parity was skipped;
- the complete 90-second cases were not run.

Important cross-language qualifications:

- NumPy and MATLAB do not naturally generate identical random sequences.
  Fixed exported signals are supported for strict comparisons.
- SciPy `firwin` matches the case's tap count, Hamming design, scaling, and
  frequency interpretation, but MATLAB-exported coefficients should be used
  when bit-level input parity is required.
- The compensation-ID recurrence follows the repository's
  `dsp.FilteredXLMSFilter` configuration and final `-flip`; exact coefficient
  parity with that System object still requires the MATLAB fixture/runtime.
- Case 2 allocates Python arrays to the actual 14-second reference length.
  MATLAB passes the full recording length to constructors but loops only over
  the shorter reference; this changes unused allocation, not recursive
  samples.
- Long recursive runs are compared by residual/weight tolerances and final
  metrics, not assumed bit-identical across runtimes.

The center-attraction term is an implemented algorithmic constraint. It is not
presented as a universal stability proof for every plant and communication
condition.
