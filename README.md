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
% Each simulation starts from a clean workspace. The search path persists
% across `clear`, but re-adding it keeps each block independently runnable.

% The simulation below will took around 15 min. t_final=1e-15 s.
clear; addpath('functions');
params_set_name='1s0_1fs';
h_atom

% The simulation below tooks around 50 minutes. t_final=1e-12 s.
clear; addpath('functions');
params_set_name='1s0_1ps';
h_atom

% The simulation below tooks around 6 hours. t_final=1e-11 s
clear; addpath('functions');
params_set_name='2p0_10ps';
h_atom

% The simulation below tooks around 6 hours. t_final=1e-11 s
clear; addpath('functions');
params_set_name='2s0_10ps';
h_atom

% Optional
% The simulation below tooks around 50 min. t_final=1e-12 s
clear; addpath('functions');
params_set_name='2p_m1_1ps';
h_atom

% Matching comparison case with m=0
clear; addpath('functions');
params_set_name='2p_m0_1ps';
h_atom

% Same density as m=1 but opposite circulation
clear; addpath('functions');
params_set_name='2p_mn1_1ps';
h_atom
```

# Running the cutoff-velocity sweep

The cutoff sweep scans the velocity cap `v_max/c` for the nodal states
`(n,l,m)=(2,1,0)` and `(2,0,0)`.  The scan runs use `M=100` trajectories and
`traj_points=1`; the integration step `dt` is given explicitly in the set name
(here `10zs` = `1e-20 s`, giving `n_steps=1e9`).  They store the energy-error
time series, nodal crossing counts, occupation fractions, and cutoff
engagement fractions in the output `.mat` file.

First run a short smoke test:

```Matlab
clear; addpath('functions');
params_set_name='test_scan_vmax_1p0c_dt_10zs';
h_atom
```

Then run the `(2,1,0)` sweep one point at a time:

```Matlab
clear; addpath('functions');
params_set_name='2p0_scan_vmax_0p01c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2p0_scan_vmax_0p03c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2p0_scan_vmax_0p1c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2p0_scan_vmax_0p3c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2p0_scan_vmax_1p0c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2p0_scan_vmax_2p0c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2p0_scan_vmax_3p0c_dt_10zs';
h_atom
```

Run the `(2,0,0)` sweep one point at a time:

```Matlab
clear; addpath('functions');
params_set_name='2s0_scan_vmax_0p01c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2s0_scan_vmax_0p03c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2s0_scan_vmax_0p1c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2s0_scan_vmax_0p3c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2s0_scan_vmax_1p0c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2s0_scan_vmax_2p0c_dt_10zs';
h_atom

clear; addpath('functions');
params_set_name='2s0_scan_vmax_3p0c_dt_10zs';
h_atom
```

For parallel runs, start separate MATLAB sessions and set one
`params_set_name` per session, for example:

```Matlab
clear; addpath('functions');
params_set_name='2p0_scan_vmax_1p0c_dt_10zs';
h_atom
```

The scan presets keep `traj_points=1` because the sweep is a statistical
energy-error test, not a trajectory-visualisation run.

Multi-trajectory runs (`M > 1`, i.e. all `*_scan_*` presets) use an
independent-seed design: the RNG seed is set by `random_seed` (default `7`)
and is appended to the run folder name as the last token, e.g.
`data/2s0_scan_vmax_1p0c_dt_10zs_seed_101/`. The scan folder name is rebuilt
from the effective `v_max/c` and `dt`, so it always carries both tags. Each
seed is stored in its own folder, so seed replicates never overwrite one
another. (Single-trajectory runs, `M == 1`, always use the fixed seed `7` and
are left untagged.)

For example, to test the `(2,0,0)` state at `v_max/c = 0.1, 1.0, 3.0` with
seed 101, run the following commands in MATLAB:

```Matlab
clear; addpath('functions');
random_seed = 101;   % set after each clear so the seed survives
params_set_name='2s0_scan_vmax_0p1c_dt_10zs';
h_atom

clear; addpath('functions');
random_seed = 101;
params_set_name='2s0_scan_vmax_1p0c_dt_10zs';
h_atom

clear; addpath('functions');
random_seed = 101;
params_set_name='2s0_scan_vmax_3p0c_dt_10zs';
h_atom
```

These three runs write to:

```text
data/2s0_scan_vmax_0p1c_dt_10zs_seed_101/
data/2s0_scan_vmax_1p0c_dt_10zs_seed_101/
data/2s0_scan_vmax_3p0c_dt_10zs_seed_101/
```

Repeat the same three commands with another `random_seed` if additional
independent replicates are needed.

# Running the time-step (dt) sweep

To check that the results are converged with respect to the integration time
step, sweep the `dt` tag of the unified scan sets at the production cutoff
`v_max/c = 1.0`.  Like the cutoff sweep these use `M=100` trajectories and
`traj_points=1` (a statistical energy test, not a trajectory-visualisation
run).  The **total simulated time is held fixed** (equal physical time across
all sweep points): `n_steps` is rescaled automatically from `dt`, so every
point covers the same physical duration.

Both `v_max/c` and `dt` are encoded in the set name, and `dt` is required.  `dt`
is given in zeptoseconds (`1 zs = 1e-21 s`), with `p` for the decimal point.
For example `_dt_100zs` → `1e-19 s`, `_dt_1zs` → `1e-21 s`, `_dt_0p5zs` →
`5e-22 s`.  For `2p0`/`2s0` the fixed total time is `1e-11 s`, so `dt=10 zs`
reproduces the production `n_steps=1e9`.

First run a short smoke test (total time `1e-13 s`):

```Matlab
clear; addpath('functions');
params_set_name='test_scan_vmax_1p0c_dt_10zs';
h_atom
```

Then run the `(2,1,0)` dt sweep one point at a time — about six logarithmically
spaced steps around the default `10 zs` (`random_seed` is set after each clear
so it survives, tagging the multi-trajectory runs):

```Matlab
clear; addpath('functions'); random_seed = 101;
params_set_name='2p0_scan_vmax_1p0c_dt_100zs';   % coarsest, n_steps=1e8
h_atom

clear; addpath('functions'); random_seed = 101;
params_set_name='2p0_scan_vmax_1p0c_dt_50zs';    % n_steps=2e8
h_atom

clear; addpath('functions'); random_seed = 101;
params_set_name='2p0_scan_vmax_1p0c_dt_10zs';    % baseline, n_steps=1e9
h_atom

clear; addpath('functions'); random_seed = 101;
params_set_name='2p0_scan_vmax_1p0c_dt_5zs';     % n_steps=2e9
h_atom

clear; addpath('functions'); random_seed = 101;
params_set_name='2p0_scan_vmax_1p0c_dt_2zs';     % n_steps=5e9
h_atom

clear; addpath('functions'); random_seed = 101;
params_set_name='2p0_scan_vmax_1p0c_dt_1zs';     % finest, n_steps=1e10
h_atom
```

The `(2,0,0)` sweep is the same with the `2s0_` prefix:

```Matlab
clear; addpath('functions'); random_seed = 101;
params_set_name='2s0_scan_vmax_1p0c_dt_100zs';
h_atom
% ... 50zs, 10zs, 5zs, 2zs, 1zs (each preceded by: clear; addpath('functions'); random_seed = 101;) ...
clear; addpath('functions'); random_seed = 101;
params_set_name='2s0_scan_vmax_1p0c_dt_1zs';
h_atom
```

These runs write to their own folders, e.g. (with `random_seed = 101`):

```text
data/2p0_scan_vmax_1p0c_dt_100zs_seed_101/
data/2p0_scan_vmax_1p0c_dt_50zs_seed_101/
data/2p0_scan_vmax_1p0c_dt_10zs_seed_101/
data/2p0_scan_vmax_1p0c_dt_5zs_seed_101/
data/2p0_scan_vmax_1p0c_dt_2zs_seed_101/
data/2p0_scan_vmax_1p0c_dt_1zs_seed_101/
```

**Runtime note:** because the total simulated time is held fixed, the step
count (and therefore the wall-clock time) scales as `1/dt`.  Halving `dt`
doubles `n_steps`; the `1 zs` point runs ~10x longer than the `10 zs`
baseline.  For parallel runs, start separate MATLAB sessions and set one
`params_set_name` per session.

To add a sweep point, just use a new `dt` (or `v_max/c`) tag in the set name —
no edit to `parameters.m` is required; both values are parsed from the name.

## Quick dt / v_max override

`dt_override` and `vmax_override` change the integration step or the velocity
cutoff at run time without defining a new set name.  Set either (or both) in the
workspace before calling `h_atom`:

- `dt_override` (seconds) replaces `dt` and rescales `n_steps` to keep the
  physical duration fixed.
- `vmax_override` (in units of `c`) replaces `v_max/c`.

For a **scan set** the output folder name is rebuilt from the *effective*
(post-override) `v_max/c` and `dt`, so it always reflects what actually ran and
never overwrites the un-overridden run:

```Matlab
clear; addpath('functions');
params_set_name='2s0_scan_vmax_1p0c_dt_10zs';
dt_override = 5e-21;     % -> dt tag 5zs
vmax_override = 3.0;     % -> vmax tag 3p0c
h_atom
% writes data/2s0_scan_vmax_3p0c_dt_5zs_seed_7/
% (the next run's leading `clear` removes dt_override / vmax_override)
```

For a **non-scan preset** (e.g. `2p0_10ps`) the override appends a `_dt_<d>zs`
(and/or `_vmax_<v>c`) tag to the preset name so the baseline run is not
overwritten.  Note this keeps the preset's own `M` and `traj_points` (for
`2p0_10ps` that is `M=1` and a large `traj_points`, i.e. a single-trajectory
visualisation run, not the `M=100` statistical sweep above):

```Matlab
clear; addpath('functions');
params_set_name='2p0_10ps';
dt_override = 5e-21;
h_atom
% writes data/2p0_10ps_dt_5zs/
% (the next run's leading `clear` restores the preset's built-in dt)
```

# Post processing
The simulation will creata 'data' directory for every parameters set. 
The post processing tooks around 5 min.
```Matlab
clear; addpath('functions');
h_atom_post_processing
```

For cutoff sweeps, the post-processing can be run directly after all sweep
points for a state have finished:

```Matlab
clear; addpath('functions');
analyse_cutoff_scan('2p0', 'figures');
analyse_cutoff_scan('2s0', 'figures');
```

To analyse a single seed only, pass the seed-specific run filter after all
three runs for that seed finish:

```Matlab
clear; addpath('functions');
analyse_cutoff_scan('2s0', 'figures/seed_101', ...
                    'data', 'seed_101');
```

This reads the folders matching `data/2s0_scan_vmax_*_seed_101` and writes:

```text
figures/seed_101/cutoff_scan_kinetic_2s0_seed_101.pdf
figures/seed_101/cutoff_scan_crossings_2s0_seed_101.pdf
```

To pool every seed available under `data/` on the same plots, use `'all'`
(the default filter loads all seeds as well):

```Matlab
clear; addpath('functions');
analyse_cutoff_scan('2s0', 'figures/all_seeds', 'data', 'all');
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
