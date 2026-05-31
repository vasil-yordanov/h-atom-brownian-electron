# Description
This repository contains the source code for simulating the hydrogen atom using 
Stochastic Mechanics to model the electron’s behavior through Brownian motion.

For more details on the methodology, see the related paper on arXiv: 
[Hydrogen Atom Simulation through Stochastic Mechanics](https://arxiv.org/abs/2412.19918v1).

# Machine Specifications

The simulation was tested on the following system:

- **Processor**: Intel(R) Core(TM) i5-10400F CPU @ 2.90GHz (6 Cores, 12 Threads)
- **RAM**: 16.0 GB

The timing results provided below are specific to this machine configuration
and may vary on systems with different specifications.

# Add functions directory to the path
Assuming you are in the home directory of the project execute:
```Matlab
addpath('functions');
```

# Running the simulation

```Matlab
% The simulation below will took around 15 min. t_final=1e-15 s.
params_set_name='1s0_1fs';
h_atom

% The simulation below tooks around 50 minutes. t_final=1e-12 s.
params_set_name='1s0_1ps';
h_atom

% The simulation below tooks around 6 hours. t_final=1e-11 s
params_set_name='2p0_10ps';
h_atom

% The simulation below tooks around 6 hours. t_final=1e-11 s
params_set_name='2s0_10ps';
h_atom

% Optional
% The simulation below tooks around 50 min. t_final=1e-12 s
params_set_name='2p_m1_1ps';
h_atom

% Matching comparison case with m=0
params_set_name='2p_m0_1ps';
h_atom

% Same density as m=1 but opposite circulation
params_set_name='2p_mn1_1ps';
h_atom
```

# Running the cutoff-velocity sweep

The cutoff sweep scans the velocity cap `v_max/c` for the nodal states
`(n,l,m)=(2,1,0)` and `(2,0,0)`.  The scan runs use `M=5` trajectories,
`dt=1e-20 s`, and `n_steps=1e9` by default.  They store the energy-error
time series, nodal crossing counts, occupation fractions, and cutoff
engagement fractions in the output `.mat` file.

First run a short smoke test:

```Matlab
addpath('functions');
params_set_name='test_scan_vmax_1p0c';
h_atom
```

Then run the `(2,1,0)` sweep one point at a time:

```Matlab
addpath('functions');
params_set_name='2p0_scan_vmax_0p01c';
h_atom

params_set_name='2p0_scan_vmax_0p03c';
h_atom

params_set_name='2p0_scan_vmax_0p1c';
h_atom

params_set_name='2p0_scan_vmax_0p3c';
h_atom

params_set_name='2p0_scan_vmax_1p0c';
h_atom

params_set_name='2p0_scan_vmax_2p0c';
h_atom

params_set_name='2p0_scan_vmax_3p0c';
h_atom
```

Run the `(2,0,0)` sweep one point at a time:

```Matlab
addpath('functions');
params_set_name='2s0_scan_vmax_0p01c';
h_atom

params_set_name='2s0_scan_vmax_0p03c';
h_atom

params_set_name='2s0_scan_vmax_0p1c';
h_atom

params_set_name='2s0_scan_vmax_0p3c';
h_atom

params_set_name='2s0_scan_vmax_1p0c';
h_atom

params_set_name='2s0_scan_vmax_2p0c';
h_atom

params_set_name='2s0_scan_vmax_3p0c';
h_atom
```

For parallel runs, start separate MATLAB sessions and set one
`params_set_name` per session, for example:

```Matlab
addpath('functions');
params_set_name='2p0_scan_vmax_1p0c';
h_atom
```

The scan presets keep `traj_points=1` because the sweep is a statistical
energy-error test, not a trajectory-visualisation run.

By default, the sweep uses a common-random-number design: each cutoff value
starts from the same seed, so differences between cutoff values are less
contaminated by unrelated Brownian histories. To run an independent-randomness
robustness check without overwriting the common sweep, set
`independent_random_design=true` and choose a `random_seed`. These runs are
written as separate run folders, for example
`data/2s0_scan_independent_seed_101_vmax_1p0c/`.

For example, to test the `(2,0,0)` state at `v_max/c = 0.1, 1.0, 3.0` with
one independent seed, run the following commands in MATLAB:

```Matlab
addpath('functions');
independent_random_design = true;
random_seed = 101;

params_set_name='2s0_scan_vmax_0p1c';
h_atom

params_set_name='2s0_scan_vmax_1p0c';
h_atom

params_set_name='2s0_scan_vmax_3p0c';
h_atom
```

These three runs write to:

```text
data/2s0_scan_independent_seed_101_vmax_0p1c/
data/2s0_scan_independent_seed_101_vmax_1p0c/
data/2s0_scan_independent_seed_101_vmax_3p0c/
```

Repeat the same three commands with another `random_seed` if additional
independent replicates are needed.

# Running the time-step (dt) sweep

To check that the results are converged with respect to the integration time
step, run the dedicated `*_scan_dt_*` parameter sets.  Like the cutoff sweep
they use `M=100` trajectories and `traj_points=1` (a statistical energy test,
not a trajectory-visualisation run), with `v_max/c` held at the production
value `1.0`.  The **total simulated time `n_steps*dt` is held fixed** (equal
physical time across all sweep points): `n_steps` is rescaled automatically from
the encoded `dt`, so every point covers the same physical duration.

The `dt` value is encoded in the set name: the mantissa uses `p` for the decimal
point, `e` separates the exponent, and `m`/`p` after `e` is the exponent sign.
For example `1em20` → `1e-20 s`, `5em21` → `5e-21 s`, `1p5em20` → `1.5e-20 s`.

The sweep is defined for the nodal states `(2,1,0)` (`2p0_scan_dt_*`) and
`(2,0,0)` (`2s0_scan_dt_*`).  For `2p0`/`2s0` the fixed total time is `1e-11 s`,
so `dt=1e-20` reproduces the production `n_steps=1e9`.

First run a short smoke test (total time `1e-13 s`):

```Matlab
addpath('functions');
params_set_name='test_scan_dt_1em20';
h_atom
```

Then run the `(2,1,0)` dt sweep one point at a time — about six logarithmically
spaced steps around the default `1e-20 s`:

```Matlab
addpath('functions');
params_set_name='2p0_scan_dt_1em19';   % coarsest, n_steps=1e8
h_atom

params_set_name='2p0_scan_dt_5em20';   % n_steps=2e8
h_atom

params_set_name='2p0_scan_dt_1em20';   % baseline, n_steps=1e9
h_atom

params_set_name='2p0_scan_dt_5em21';   % n_steps=2e9
h_atom

params_set_name='2p0_scan_dt_2em21';   % n_steps=5e9
h_atom

params_set_name='2p0_scan_dt_1em21';   % finest, n_steps=1e10
h_atom
```

The `(2,0,0)` sweep is the same with the `2s0_` prefix:

```Matlab
addpath('functions');
params_set_name='2s0_scan_dt_1em19';
h_atom
% ... 5em20, 1em20, 5em21, 2em21, 1em21 ...
params_set_name='2s0_scan_dt_1em21';
h_atom
```

These runs write to their own folders, e.g.:

```text
data/2p0_scan_dt_1em19/
data/2p0_scan_dt_5em20/
data/2p0_scan_dt_1em20/
data/2p0_scan_dt_5em21/
data/2p0_scan_dt_2em21/
data/2p0_scan_dt_1em21/
```

**Runtime note:** because the total simulated time is held fixed, the step
count (and therefore the wall-clock time) scales as `1/dt`.  Halving `dt`
doubles `n_steps`; the `1e-21 s` point runs ~10x longer than the `1e-20 s`
baseline.  For parallel runs, start separate MATLAB sessions and set one
`params_set_name` per session.

To add a sweep point, just use a new encoded `dt` in the set name — no edit to
`parameters.m` is required; the value is parsed from the name.

## Quick dt override on any preset

For a one-off check on a preset that has no `scan_dt` set, set `dt_override` in
the workspace before calling `h_atom`.  It replaces `dt` and rescales `n_steps`
the same way (equal physical time), and writes to a folder tagged with the `dt`
value in zeptoseconds (`1 zs = 1e-21 s`, `p` for the decimal point), e.g.
`dt_override = 5e-21` -> `data/2p0_10ps_dt_5zs/`, so the baseline run is not
overwritten.  Note this keeps the preset's own `M` and `traj_points` (for
`2p0_10ps` that is `M=1` and a large `traj_points`, i.e. a single-trajectory
visualisation run, not the `M=100` statistical sweep above).

```Matlab
addpath('functions');
params_set_name='2p0_10ps';
dt_override = 5e-21;
h_atom
clear dt_override   % restore the preset's built-in dt afterwards
```

# Post processing
The simulation will creata 'data' directory for every parameters set. 
The post processing tooks around 5 min.
```Matlab
h_atom_post_processing
```

For cutoff sweeps, the post-processing can be run directly after all sweep
points for a state have finished:

```Matlab
addpath('functions');
analyse_cutoff_scan('2p0', 'figures');
analyse_cutoff_scan('2s0', 'figures');
```

For an independent-randomness robustness run, use the same analysis function
with the seed-specific run filter after all three independent runs finish:

```Matlab
addpath('functions');
analyse_cutoff_scan('2s0', 'figures/independent_seed_101', ...
                    'data', 'independent_seed_101');
```

This reads the independent-seed folders matching
`data/2s0_scan_independent_seed_101_vmax_*` and writes:

```text
figures/independent_seed_101/cutoff_scan_kinetic_2s0_independent_seed_101.pdf
figures/independent_seed_101/cutoff_scan_crossings_2s0_independent_seed_101.pdf
```

To overlay the common-random sweep and all independent-seed points available
under `data/` on the same plots, use:

```Matlab
addpath('functions');
analyse_cutoff_scan('2s0', 'figures/all_random_designs', 'data', 'all');
```

In this combined plot, the common-random sweep is drawn as connected curves,
while independent-seed runs are overlaid as unconnected markers. This allows
independent checks with different seeds at different cutoff values without
mixing them into one artificial sweep curve.

This writes the cutoff sensitivity and crossing-count figures to `figures/`,
including:

```text
figures/cutoff_scan_kinetic_2p0.pdf
figures/cutoff_scan_crossings_2p0.pdf
figures/cutoff_scan_kinetic_2s0.pdf
figures/cutoff_scan_crossings_2s0.pdf
```

If `h_atom_post_processing` is called with `params_set_name` set to a scan
run, it also checks whether all cutoff points for that state are present and
regenerates the scan figures when the sweep is complete.

# Results

**Note:** Due to the performance of the movie, we do not show the detailed trajectory of the particle
for later moments (see [plot_trajectory_3d.m](functions/plot_trajectory_3d.m)), the algorithm is:
```
  if ~show_all
        % Calculate the step size, ensuring it's at least 1
        every = max(floor(traj_points / 1e6), 1);
        idx = 1:every:m;
    else
        idx = 1:m;
    end
```

## Brownian motion on sphere with Bohr radius
[![Brownian motion on sphere](movies/1s0_1fs.gif)](movies/1s0_1fs.gif)

## $(n,l,m)=(1,0,0)$ state
[![1s0_state](movies/1s0_1ps.gif)](movies/1s0_1ps.gif)

##  $(n,l,m)=(2,1,0)$ state
[![2p0_state](movies/2p0_10ps.gif)](movies/2p0_10ps.gif)

## $(n,l,m)=(2,1,1)$ state
[![2p1_state](movies/2p1_1ps.gif)](movies/2p1_1ps.gif)

Recommended preset labels: `2p_m0_1ps`, `2p_m1_1ps`, `2p_mn1_1ps`. The legacy names `2p0_1ps` and `2p1_1ps` remain available as aliases.


## $(n,l,m)=(2,0,0)$ state
[![2s0_state](movies/2s0_10ps.gif)](movies/2s0_10ps.gif)
