# BLAKJac

[![Build Status]([https://github.com/Miha Fuderer/BLAKJac.jl](https://github.com/mfuderer/BLAKJac.jl)/actions/workflows/CI.yml/badge.svg?branch=main)]([https://github.com/Miha Fuderer/BLAKJac.jl](https://github.com/mfuderer/BLAKJac.jl)/actions/workflows/CI.yml?query=branch%3Amain)

## Introduction
BLAKJac is implemented in the Julia language. It relates to the MR-STAT functionality and aims to predict or minimize the noise levels in the reconstructed quantitative maps.

Two of its main functions are
* Analyzer: given a sequence of encoding values and of RF pulse angles, it predicts the resulting noise level.
* Optimizer: given a sequence of encoding values, it generates a sequence of RF pulse angles that is presumed to be optimal in terms of noise level.
The Optimizer makes use of the Analyzer.

## Inputs and outputs
### Analyzer
In essence, the analyzer has as its input
* The RF sequence: a vector of $N_{TR}$ complex values, with e.g. $N_{TR}=5\cdot 224=1120$. The magnitude of the complex values expresses the flip angle in degrees.
* The *trajectorySet*. See dedicated section below.
* A large bunch of options, in the form of a dictionary. For details, see BLAKJac_analysis and BLAKJac_defaults

There are two additional parameters:
* The resource: what type of processor is available;
* [Actually an output] Optionally, a dictionary of $I^{-1}$ matrices, where $I$ is the diagonalized Fisher Information Matrix. Such a dictionary may be useful if the optimizer is run a multiplicity of times.

As an *Output*, the Analyzer returns
* A vector of noise levels of length $p$, with $p$ being the number of properties. For the time being, $p$ is always 3.
* An estimate of the *information content* that can be expected from the scan; for the time being, a pretty immature metric.
* *b1factorRMS*: an indication on how sensitive the sequence is to deviations of B1+.

Mostly, only the first output is used, the vector of noise levels. These have a somewhat quantitative meaning: these are scaled such that if $T_1$ and $T_2$ were to be known beforehand, an ideal sequence would result in the vector [1.0, 0.0, 0.0]. Typically, the results are on the order of 3; a value of 2 is “good”, a value over 4 is generally bad and if it results in values below 1.0, then something is fishy – unless strong regularization has been applied (see “options”).

The interface, in full:

```julia
function BLAKJac_analysis!(
    resource::CPU1, 
    RFdeg::Vector{ComplexF64}, 
    trajectorySet::Vector{Vector{TrajectoryElement}}, 
    options::Dict, 
    saved_H::Dict=Dict())
```

### Optimizer
In brief, the optimizer has a trajectorySet as input and it spits out an RF sequence. And it needs a lot of options. 

The input:
* A *trajectorySet*. See dedicated section below.
* A large bunch of options, in the form of a dictionary.
* (Optionally) a random-value seed. Relevant if the option opt_initialize (see below) is set to cRandom30. By default, every run of the optimizer will always return the same pattern. If different realizations are needed, each call to the optimizer must be done with a different seed (Int).

Output
* It returns an RF sequence: a vector of $N_{TR}$ complex values, with e.g. $N_{TR}=5\cdot 224=1120$. The magnitude of the complex values expresses the flip angle in degrees.

The interface, in full:

function BLAKJac_optimize(trajectorySet, options::Dict, mersenneTwisterSeed=1)

## TrajectorySet
The *trajectorySet*. This is a vector of $N_{TR}$ vectors (in Cartesian applications, these are vectors of length 1) of *TrajectoryElement*. And a TrajectoryElement is a tuple, mostly used for $(k_y,k_z)$ – and then, in 2D-cases, the value of $k_z$ is fixed (typically, always 0).

Note that in the 2D case, $k_x$ is actually ignored, assuming that all samples of the readout share the same history – which ignores T2-decay during readout.

The reason of making the trajectorySet a vector of vectors is for non-Cartesian: then, it can be e.g. a vector of 20 spirals of 1000 samples each and “$(k_y,k_z)$” has to be interpreted as $(k_x,k_y)$.

The function TrajectorySet(kyVector, kzVector) may be useful before invoking the Analyzer: it converts an $N_{TR}$-vector of $k_y$ values and an $N_{TR}$-vector of $k_z$ values into an $N_{TR}$-vector of singleton-vectors of tuples $(k_y,k_z)$. That function has the following (tricky?) logic: if any of the provided $k_y$ or $k_z$ values is negative, then it assumes $(k_y,k_z)=(0,0)$ corresponds to zero spatial frequency. If the provided values are all nonnegative, it assumes the center of the range to correspond to zero spatial frequency. This logic is due to legacy reasons, where typically $k_y=1:224$ and $k_z=1$ would be provided, so TrajectorySet shifts the values by (112,1) in that example.

## Options
There are lots of option parameters, some of them are deprecated and some have defaults that are very rarely adapted. For that purpose, the following function is available:

BLAKJac_defaults!(trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)

It does not return anything, but the “!” indicates that options will be modified by it. Note: it only assigns uninitialized options to a default value. Options that are already filled in will be left unmodified. To set all options to default, use

recon_options = Dict() # erase all existing settings 

BLAKJac.BLAKJac_defaults!(trajectorySet, recon_options)

The *trajectorySet* is a parameter, since some defaults are set accordingly to the encoding sequence.

The list of options:

| Name | Intended type | Default | Explanation |
| :---- | :---- | :---- | :---- |
| rfFile | string | (empty) | Functionally disused. Its only purpose is to tag figures (which may be generated inside the Analyzer) with the selected string. |
| rflabel | string | (empty) | Functionally disused. Its only purpose is to tag figures (which may be generated inside the Analyzer) with the selected string. |
| nTR | Int | Number of trajectories | Number of excitation pulses |
| **TR** | Float | 0.01 | TR in seconds | 
| T1ref | Float | 0.67 | The reference T1 value - by now, mainly used for noise scaling and for "relative" noise deviations |
| T2ref | Float | 0.076 | The reference T2 value - by now, mainly used for noise scaling and for "relative" noise deviations |
| **startstate** | Float | 1.0 | Although it is formally a float, it acts like a Boolean: it indicates the starting state of the z-magnetisation; +1 for no prepulse, -1 for inversion |
| maxstate | Int | 64 |  Length of history taken into account in simulating magnetisation from sequence |
| nky | Int | Derived from range of ky values in trajectorySet | Length of range of ky encoding values. |
| nkz | Int | Derived from range of kz values in trajectorySet | Length of range of kz encoding values. So nTR/nky/nkz is number of samples per encoding. | 
| maxMeas | Int | Filled in as 4*nTR/(nky*nkz) | (R1) (see remarks below) |
| handleB1 | string | “no” | Available for testing the B1-sensitivity of sequences. Three valid values are "no", "sensitivity", "co-reconstruct". If “sensitivity” is selected, then the value of b1factorRMS is a meaningful output of the Analyzer: a value of 0 indicates that the T1 or T2 maps will not be affected by a deviation of B1; an output of e.g. -1.5 for e.g. T2 (which is typical) indicates that 1% of error in B1 will cause a 1.5% bias (underestimation) of T2. With “co-reconstruct”, B1 is just considered as a 4th parameter to be reconstructed (alongside ρ, T1 and T2) |
| considerCyclic | Bool | false | If true, it is assumed that the RF sequence repeats itself endlessly without pause. Only compatible with startstate==1 |
| useSymmetry | Bool | false | (Oops. Blooper in default setting: should have been made “true”); Assume real-ness of rho, T1 and T2 and therefor symmetry in k-space. |
| invregval | Vector{Float} | [200.0^(-2), 200.0^(-2), 200.0^(-2)] | (R2) |
| useSurrogate | Bool | false | (Relic. Don’t touch) |
| sigma_ref | Float | 0.2 | (Irrelevant unless studying the information-content criterion, yet immature). Reference normalized noise level for information-content estimation. See logbook around 2021-06-03 |
| sar_limit | Float | 40^2/0.01 | SAR limit. Initially set to an RMS level of 40 degrees at TR=10ms | 
| T1T2set | Array{(Float,Float)} | (R3) | A set of (e.g. 7) points in (T1,T2)-space around which the Analyzer will evaluate the performance. The output is the average over the performances at each point. The current default has been set to span a clinically relevant range of relaxation values. It has been set with 1.5T and 3T in mind. For very different fields (e.g. 7T), adaptation of the T1 values may be indicated. |
| sizeSteps | Vector{Int} | [10] | (R4) |
| portion_size | Int | 50 | (R4) |
| portion_step | Int | 30 | (R4) |
| opt_iterations_limit | Int | (R4) | (R4) |
| optpars | Optim.Options | (R5) | (R5) |
| opt_method |  | NelderMead | Optimizer algorithm. Do not change: most other choices require the explicit knowledge of the gradient of the criterion, or require a prohibitive number of iterations. |
| opt_expand_type | string | “spline” | Deprecated. Do not change |
| opt_keep_positive | Bool | False | If true, resulting RF angles are not allowed to be negative. |
| opt_initialize | string | “cRandom30” | Optimizer option: defines the starting pattern for the optimization process. Valid values are: "cRandom30", "ones", "ernst", "init_angle", "quadraticPhase30", "RF_shape". With cRandom30, a random pattern of angles, with standard deviation of 30 degrees, is chosen. “ones” and “ernst” initialize to a flat initialization. With the choice “RF_shape”, it is initialized by the function rfFunction (see below) |
| rfFunction | function | (undef) | A function with the parameters nTR and nky that outputs a vector of RF pulse values. | 
| opt_complex | Bool | True | Optimizer results in a complex output. (This option is actually deprecated, since the following one gives much better results; so the current default is strange: has to be changed to false) |
| opt_slow_phase | Bool | false | The phase of the pulses becomes an integral of an integral of a pattern that is subject to optimization alongside the amplitude (or, strictly speaking: the real part) of the RF pulses. |
| opt_imposed_2nd_derivative_of_phase | Vector{Float} | [] | Allows to overrule aforementioned optimized pattern. If empty, no overruling takes place. So it can be used for a “hand-drawn” 2nd derivative of the phase, provided in rad. The length of the vector must equal nTR. |
| opt_focus | String | “mean” | Valid values are: "rho","T1","T2","mean","weighted", "max". If e.g. selecting “T2”, then the sequence will optimize on the noisiness of the T2 map. Despite its default value, currently, the “max” is en vogue. This selects the criterion $\max({\sigma_{T1}}/{T_{1,ref}},{\sigma_{T2}}/{T_{2,ref}})$. |
| emphasize_low_freq | Bool | false | if switched on, then the optimization will steer more on reducing the noise factor on low spatial frequencies. | 
| opt_criterion | string | "sar-limited noise" | (R6) |
| account_SAR | Bool | false | SAR-limit has to be taken into account |
| lambda_B1 | Float | 10.0 | penalty for B1 sensitivity |
| lambda_CSF | Float | 0.0 | penalty for CSF sensitivity |
| opt_emergeCriterion | Int | 1000 | Every opt_emergeCriterion iterations, the Optimizer will display an intermediate result. |
| plottypes | Vector of string | [] | Controls (debugging) plots. Can be a collection of the following values: "first", "trajectories", "weights", "angles", "noisebars", "noiseSpectrum", "infocon" |
| optcount | Int | 0 | (For internal use – do not change) |
| stage_text | string | “” | (For internal use – do not change) |
| plotfuncs | Dict{functions} |  | (For internal use – do not change) |

Remarks: 

(R1)
: maxMeas: Maximum expected number of measurements per phase-encoding value. A ToDo item: this should be automatic, but for the time being, it has to be provided by the user. If there is no variable density, then the default value has to be set to 4 times the number of k-spaces, i.e. 4*nTR/nky in 2D sequences and to 8*nTR/(nky*nkz) if two phase-encoding directions are applied (which points to a flaw in the default, in the case of 3D).

Why 4 or 8 and not 1? One factor of 2 is needed because of assumed positivity of T1 and T2, for which the signal at $(k_y,k_z)$ and at $(-k_y,-k_z)$ has to be considered jointly. Another factor of 2 (per dimension) is needed because the entered $k_y$ or $k_z$ values are allowed to be fractional.
In the case of a fractional value of $k_y$, e.g. at $k_y=n+f$ with $f\in [0,1)$ and $n$ an integer, part of the signal is assigned to the integer location $n$ and another part to integer location $n+1$. To be precise, a fraction $\cos({\pi f}/{2})$ is assigned to $n$ and $\sin({\pi f}/{2})$ to $(n+1)$. This implies that we have some, albeit minimalistic, multi-coil capacity in our system. The cos() and sin() functions have been chosen in order to make this “splitting” noise-neutral, since it preserves the total noise power.

(R2)
: invregval: Although the actual reconstruction does not apply explicit regularization, the Analyzer assumes that it does. With the assumed default of $200^{-2}$, the Analyzer acts as if regularization would bias the reconstructed maps towards $(\rho,T_1,T_2)=(0,T_{1,ref},T_{2,ref})$, by striving for a minimal squared error in the case that the SNR is around 200. In our case, “SNR” has to be seen as the expected value of ${\sqrt{(1/N_{voxels})\sum_{voxels}(T_1(r)-T_{1,ref})^2}}/{\sigma_{T_1}}$.
The parameter can also be used to indicate that a given parameter is “known”. E.g. by setting invregval to $[10^{+6}.10^{+6},200^{-2}]$, we indicate that  and T1 are known and that all data are assumed to be used to estimate T2. 

(R3)
: The current default of the T1T2set is set to
[(0.8183, 0.0509), (0.67, 0.2066), (1.2208, 0.1253), (0.2465, 0.0461), (2.2245, 0.3082), (0.3677, 0.0461), (0.3677, 0.1253)]

(R4)
: The options sizeStep, portion_size, portion_step and opt_iterations_limit are Opimizer options. Following approach is used: in a first stage the Vector sizeStep is relevant; let it be $[n_1,n_2,n_3]$ as an example. In that first stage, the Optimizer internally runs a series of optimizations. The first one is to optimize the flip angle on $n_1$ equidistant points, the other flip angles following by cubic-spline interpolation. Subsequently, the result is cubic-spline interpolated to $n_2$ equidistant points to act as a starting point on an optimization on $n_2$ points, etc. Once this is done, the result is used as a starting point for the second stage of optimizations (deprecated by now): the Optimizer will then optimize *individual* RF pulses. Since the optimization process cannot reasonably handle a large number (over 50 or so) of parameters simultaneously, it only addresses a portion of size portion_size at a time, taking all other RF pulses as fixed. Once such a portion is considered optimized, the portion advances by portion_step pulses. In some cases, this may result in sequences with pulses around 180 degrees, but not if a sar_limit has been imposed to e.g. 40 degrees RMS.
The total number of internal runs should therefore be

 length(sizeStep)+ceil(nTR/portion_step). 

This can be reduced by the option parameter opt_iterations_limit. When omitting the second stage (which is advised), the default opt_iterations_limit=length(sizeStep) is appropriate.

(R5)
: Internal Optimizer options. Default is 
Optim.Options(time_limit = 3000.0, iterations = 100000, f_tol=1.0e-2, g_tol = 1.0e-3)
The first two options are useful to stop the prevent “the optimizer being stuck, either by non-convergence or due to a parameter choice that would require ages to converge. 
On f_tol and g_tol, the latest insights are that these should be set stricter in order to avoid grossly-suboptimal false minima. A value of around $10^{-8}$ should be default.

(R6)
: Valid values are:
* "noise level" (not accounting for SAR or high flip angles), 
* "low-flip corrected noise" (penalty on the maximum flip angle in the sequence, on the assumption that, on a fixed TR, it will thereby reduce the sampling-time and consequently increase the noise by a certain factor), 
* "information content", "information content penalized" (immature – do not use), 
* (Legacy )"sar-limited noise" (a heavy penalty term is added if the root-mean-square of the RF pulses exceeds a limit). Deprecated. Instead, combine "noise level" with account_SAR=true




















