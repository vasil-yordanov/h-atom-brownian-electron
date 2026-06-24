# Description
This repository contains the source code for simulating the hydrogen atom using
Stochastic Mechanics to model the electron’s behavior through Brownian motion.

For more details on the methodology, see the related paper on arXiv:
[Hydrogen Atom Simulation through Stochastic Mechanics](https://arxiv.org/abs/2412.19918v1).

# Results

**Note:** Due to the performance of the movie, we do not show the detailed
trajectory of the particle for later moments (see
[plot_trajectory_3d.m](functions/plot_trajectory_3d.m)).

## Brownian motion on sphere with Bohr radius
[![Brownian motion on sphere](movies/1s0_1fs.gif)](movies/1s0_1fs.gif)

## $(n,l,m)=(1,0,0)$ state
[![1s0_state](movies/1s0_1ps.gif)](movies/1s0_1ps.gif)

The electron's initial position is intentionally chosen far from the nucleus, at
spherical coordinates $(10\,a_0,\ \pi/2,\ 0)$, to demonstrate the drift toward
the nucleus.

##  $(n,l,m)=(2,1,0)$ state
[![2p0_state](movies/2p0_10ps.gif)](movies/2p0_10ps.gif)

## $(n,l,m)=(2,1,1)$ state
[![2p1_state](movies/2p1_1ps.gif)](movies/2p1_1ps.gif)

## $(n,l,m)=(2,0,0)$ state
[![2s0_state](movies/2s0_10ps.gif)](movies/2s0_10ps.gif)

Recommended preset labels: `2p_m0_1ps`, `2p_m1_1ps`, `2p_mn1_1ps`. The legacy names `2p0_1ps` and `2p1_1ps` remain available as aliases.

## Phase-driven azimuthal circulation ($m = \pm 1, 0$)

The $m=+1$ and $m=-1$ states circulate in opposite directions while $m=0$ does
not. Looking straight down the $+z$ axis, the electron (dark dot, fading comet
tail) sweeps counter-clockwise for $m=+1$, clockwise for $m=-1$, and shows no
net sense of rotation for $m=0$:

<a href="movies/2p_pm1_azimuthal_current_3panel.gif"><img src="movies/2p_pm1_azimuthal_current_3panel.gif" width="70%"></a>

How to reproduce everything is below: first
[run the simulations](#running-the-simulations), then
[generate the README figures](#generating-the-readme-figures) (the movies above)
and the [manuscript figures](#generating-the-manuscript-figures).

# Machine Specifications

The simulation was run on a **MacBook Pro** (`MacBookPro18,1`), **Apple M1 Pro**
(10 cores: 8 performance + 2 efficiency), **16 GB** RAM, macOS 15.3.1. The run
times quoted in this README are specific to this machine and will vary on other
hardware.

# Running the simulations

Every simulation and figure is produced by a `matlab -batch` one-liner. On macOS
the binary may be at `/Applications/MATLAB_R2024b.app/bin/matlab`, so use your
own path. Each run writes a `data/<run>/` folder (its `.mat` data and, for the
production runs, an `.mp4` dashboard movie).

## Production runs (Figs 1-10 and the movies)

Single-trajectory runs (`M=1`); the seed is fixed to `7` automatically, so no
`random_seed` is needed.

```bash
matlab -batch 'addpath("functions"); params_set_name="1s0_1fs";   h_atom'   # ~16 min
matlab -batch 'addpath("functions"); params_set_name="1s0_1ps";   h_atom'   # ~30 min
matlab -batch 'addpath("functions"); params_set_name="2p0_10ps";  h_atom'   # ~2.9 h
matlab -batch 'addpath("functions"); params_set_name="2s0_10ps";  h_atom'   # ~2.9 h
matlab -batch 'addpath("functions"); params_set_name="2p_m1_1ps"; h_atom'   # ~27 min
matlab -batch 'addpath("functions"); params_set_name="2p_m0_1ps"; h_atom'   # ~33 min
matlab -batch 'addpath("functions"); params_set_name="2p_mn1_1ps"; h_atom'  # ~37 min
```

## Cutoff scan, (2,0,0) state (Fig. 11)

Statistical runs with `M=100` trajectories and `traj_points=1`, fixed total time
`T = 2e-11 s` (20 ps); `n_steps = T/dt` is set automatically. The grid is the
eight cutoffs `v_max/c = 0.01, 0.03, 0.1, 0.3, 1, 2, 3, 5` and the four steps
`dt = 5, 10, 50, 100 zs`, so 32 runs. Each line below carries the exact
`random_seed` used in the published figure.

```bash
# dt = 5 zs  (~28 h each)
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p01c_dt_5zs"; random_seed=7379; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p03c_dt_5zs"; random_seed=2757; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p1c_dt_5zs";  random_seed=1115; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p3c_dt_5zs";  random_seed=1213; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_1p0c_dt_5zs";  random_seed=1008; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_2p0c_dt_5zs";  random_seed=1005; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_3p0c_dt_5zs";  random_seed=1007; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_5p0c_dt_5zs";  random_seed=1009; h_atom'   # ~28 h

# dt = 10 zs  (~14 h each)
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p01c_dt_10zs"; random_seed=82379; h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p03c_dt_10zs"; random_seed=41757; h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p1c_dt_10zs";  random_seed=3227;  h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p3c_dt_10zs";  random_seed=3999;  h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_1p0c_dt_10zs";  random_seed=3001;  h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_2p0c_dt_10zs";  random_seed=3005;  h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_3p0c_dt_10zs";  random_seed=3007;  h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_5p0c_dt_10zs";  random_seed=3009;  h_atom'   # ~14 h

# dt = 50 zs  (~3 h each)
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p01c_dt_50zs"; random_seed=7913; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p03c_dt_50zs"; random_seed=7815; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p1c_dt_50zs";  random_seed=7773; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p3c_dt_50zs";  random_seed=7313; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_1p0c_dt_50zs";  random_seed=5213; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_2p0c_dt_50zs";  random_seed=5315; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_3p0c_dt_50zs";  random_seed=5513; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_5p0c_dt_50zs";  random_seed=5713; h_atom'   # ~3 h

# dt = 100 zs  (~1.4 h each)
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p01c_dt_100zs"; random_seed=8113;  h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p03c_dt_100zs"; random_seed=8217;  h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p1c_dt_100zs";  random_seed=8111;  h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_0p3c_dt_100zs";  random_seed=8819;  h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_1p0c_dt_100zs";  random_seed=98223; h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_2p0c_dt_100zs";  random_seed=98335; h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_3p0c_dt_100zs";  random_seed=99121; h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2s0_scan_vmax_5p0c_dt_100zs";  random_seed=99819; h_atom'   # ~1.4 h
```

## Cutoff scan, (2,1,0) state (Fig. 11 companion)

Same grid and parameters with the `2p0_` prefix.

```bash
# dt = 5 zs  (~28 h each)
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p01c_dt_5zs"; random_seed=9359;  h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p03c_dt_5zs"; random_seed=2359;  h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p1c_dt_5zs";  random_seed=67870; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p3c_dt_5zs";  random_seed=8432;  h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_1p0c_dt_5zs";  random_seed=34117; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_2p0c_dt_5zs";  random_seed=26418; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_3p0c_dt_5zs";  random_seed=87576; h_atom'   # ~28 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_5p0c_dt_5zs";  random_seed=23870; h_atom'   # ~28 h

# dt = 10 zs  (~14 h each)
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p01c_dt_10zs"; random_seed=83356; h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p03c_dt_10zs"; random_seed=37576; h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p1c_dt_10zs";  random_seed=31418; h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p3c_dt_10zs";  random_seed=31117; h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_1p0c_dt_10zs";  random_seed=1118;  h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_2p0c_dt_10zs";  random_seed=2109;  h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_3p0c_dt_10zs";  random_seed=2129;  h_atom'   # ~14 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_5p0c_dt_10zs";  random_seed=2669;  h_atom'   # ~14 h

# dt = 50 zs  (~3 h each)
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p01c_dt_50zs"; random_seed=6519; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p03c_dt_50zs"; random_seed=6577; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p1c_dt_50zs";  random_seed=6588; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p3c_dt_50zs";  random_seed=2212; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_1p0c_dt_50zs";  random_seed=7718; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_2p0c_dt_50zs";  random_seed=7878; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_3p0c_dt_50zs";  random_seed=7919; h_atom'   # ~3 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_5p0c_dt_50zs";  random_seed=9879; h_atom'   # ~3 h

# dt = 100 zs  (~1.4 h each)
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p01c_dt_100zs"; random_seed=88576; h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p03c_dt_100zs"; random_seed=87576; h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p1c_dt_100zs";  random_seed=39576; h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_0p3c_dt_100zs";  random_seed=73870; h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_1p0c_dt_100zs";  random_seed=39879; h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_2p0c_dt_100zs";  random_seed=39870; h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_3p0c_dt_100zs";  random_seed=79879; h_atom'   # ~1.4 h
matlab -batch 'addpath("functions"); params_set_name="2p0_scan_vmax_5p0c_dt_100zs";  random_seed=76879; h_atom'   # ~1.4 h
```

Each run writes to its own folder `data/<set>_seed_<seed>/`. Run points in
parallel by launching one command per MATLAB session.

## Cutoff sweep run time

The total simulated time is fixed at `T = 2e-11 s` (20 ps), so the number of
steps is `n_steps = T/dt` and the run time scales as `1/dt`, almost independent
of `v_max/c`. The measured `dt = 5 zs` runs took about 25.8 to 29.6 h (around
28.5 h on average, about 2.57e-5 s/step at `n_steps = 4e9`). Running 8 runs at a
time (one MATLAB session per core), the eight `v_max/c` points of a given `dt`
finish in about one per-run time:

| dt | n_steps | per run | one batch of 8 (parallel) |
|------|-----------|---------|---------------------------|
| 5 zs   | 4e9 | ~28.5 h | ~28.5 h |
| 10 zs  | 2e9 | ~14.3 h | ~14.3 h |
| 50 zs  | 4e8 | ~2.9 h  | ~2.9 h  |
| 100 zs | 2e8 | ~1.4 h  | ~1.4 h  |

Running the four dt batches one after another, one state (32 runs) takes about
28.5 + 14.3 + 2.9 + 1.4 = about 47 h (~2 days), and both states (64 runs) about
94 h (~3.9 days). The estimates are good to about +/-15%.

# Generating the README figures

These are the movies and image shown in [Results](#results). Run the relevant
[simulations](#running-the-simulations) first.

## State movies

Every production run writes its 6-panel dashboard as an MP4 to
`data/<run>/<run>.mp4` (this is on by default; set `make_live_plots = false` in
`h_atom.m` to skip it). Convert each MP4 to a GIF with
[ffmpeg](https://ffmpeg.org) (`brew install ffmpeg`). The dashboard is rendered
at 1200 px wide so the six panels stay legible; `dither=none` keeps the flat
background clean. The four long movies are sampled to 5 fps (to keep the GIF
size reasonable at this resolution); the short (2,1,1) run keeps its native rate:

```bash
export PAL="split[s0][s1];[s0]palettegen=max_colors=128:stats_mode=diff[p];[s1][p]paletteuse=dither=none:diff_mode=rectangle"
ffmpeg -y -i data/1s0_1fs/1s0_1fs.mp4     -vf "fps=5,scale=1200:-1:flags=lanczos,$PAL" movies/1s0_1fs.gif
ffmpeg -y -i data/1s0_1ps/1s0_1ps.mp4     -vf "fps=5,scale=1200:-1:flags=lanczos,$PAL" movies/1s0_1ps.gif
ffmpeg -y -i data/2p0_10ps/2p0_10ps.mp4   -vf "fps=5,scale=1200:-1:flags=lanczos,$PAL" movies/2p0_10ps.gif
ffmpeg -y -i data/2p_m1_1ps/2p_m1_1ps.mp4 -vf "fps=5,scale=1200:-1:flags=lanczos,$PAL" movies/2p1_1ps.gif
ffmpeg -y -i data/2s0_10ps/2s0_10ps.mp4   -vf "fps=5,scale=1200:-1:flags=lanczos,$PAL" movies/2s0_10ps.gif
```

The four long movies (1 fs / 1 ps / 10 ps simulated) play over 100 s at 5 fps.
The (2,1,1) movie is short and keeps its native frame rate. The run name and the
GIF name differ only for that state:

| State | Run (`params_set_name`) | MP4 written | GIF |
|-------|-------------------------|-------------|-----|
| Brownian sphere | `1s0_1fs`   | `data/1s0_1fs/1s0_1fs.mp4`     | `movies/1s0_1fs.gif`  |
| (1,0,0) | `1s0_1ps`   | `data/1s0_1ps/1s0_1ps.mp4`     | `movies/1s0_1ps.gif`  |
| (2,1,0) | `2p0_10ps`  | `data/2p0_10ps/2p0_10ps.mp4`   | `movies/2p0_10ps.gif` |
| (2,1,1) | `2p_m1_1ps` | `data/2p_m1_1ps/2p_m1_1ps.mp4` | `movies/2p1_1ps.gif`  |
| (2,0,0) | `2s0_10ps`  | `data/2s0_10ps/2s0_10ps.mp4`   | `movies/2s0_10ps.gif` |

### Rebuilding a movie from saved data (no re-simulation)

`make_dashboard_movie` regenerates a run's dashboard `.mp4` from the data the
simulation already saved, so you can restyle the panels and rebuild in minutes
instead of re-running the (hours-long) simulation.

First rebuild the dashboard MP4s from the saved data (MATLAB):

```bash
matlab -batch 'make_dashboard_movie("1s0_1fs")'
matlab -batch 'make_dashboard_movie("1s0_1ps")'
matlab -batch 'make_dashboard_movie("2p0_10ps")'
matlab -batch 'make_dashboard_movie("2s0_10ps")'
matlab -batch 'make_dashboard_movie("2p_m1_1ps")' 
```

Then convert the MP4s to GIFs (ffmpeg):

```bash
export PAL="split[s0][s1];[s0]palettegen=max_colors=128:stats_mode=diff[p];[s1][p]paletteuse=dither=none:diff_mode=rectangle"
ffmpeg -y -i data/1s0_1fs/1s0_1fs.mp4     -vf "$PAL" movies/1s0_1fs.gif
ffmpeg -y -i data/1s0_1ps/1s0_1ps.mp4     -vf "$PAL" movies/1s0_1ps.gif
ffmpeg -y -i data/2p0_10ps/2p0_10ps.mp4   -vf "$PAL" movies/2p0_10ps.gif
ffmpeg -y -i data/2s0_10ps/2s0_10ps.mp4   -vf "$PAL" movies/2s0_10ps.gif
ffmpeg -y -i data/2p_m1_1ps/2p_m1_1ps.mp4 -vf "$PAL" movies/2p1_1ps.gif
```

It writes the same `data/<run>/<run>.mp4` path as the live run. It is a replica
of the live rendering loop in `h_atom.m`: it redraws each frame from the exact
per-frame data the simulation stored (`hist_counts_*_traj`, the energy and
distribution-deviation time series, and the trajectory), with the same plotting
functions. The default arguments downsize to the GIF target (1200 px, 5 fps);
pass `2100, 1` for the full-resolution, full-frame-rate exact replica of the
live MP4 (`make_dashboard_movie("1s0_1fs", 2100, 1)`). The short `2p_m1_1ps` run
only has 80 frames, so it uses `frame_stride = 1` to keep them all.

**Both ways give the same movie.** The simulation stores each frame's dashboard
state, and `make_dashboard_movie` feeds exactly those snapshots back into the
same plot functions — so with `make_dashboard_movie("<run>", 2100, 1)` the result
is identical to the movie the simulation wrote live; it does not matter which way
you generate it. To confirm on your machine:

```bash
# (a) live movie written by the simulation
matlab -batch 'addpath("functions"); params_set_name="1s0_1fs"; h_atom'
cp data/1s0_1fs/1s0_1fs.mp4 /tmp/live.mp4

# (b) movie rebuilt from the saved data, exact-replica settings
matlab -batch 'make_dashboard_movie("1s0_1fs", 2100, 1)'
cp data/1s0_1fs/1s0_1fs.mp4 /tmp/rebuilt.mp4

# (c) compare frame 500 -- the max pixel difference should be ~0
ffmpeg -y -i /tmp/live.mp4    -vf "select=eq(n\,500)" -frames:v 1 /tmp/live_500.png
ffmpeg -y -i /tmp/rebuilt.mp4 -vf "select=eq(n\,500)" -frames:v 1 /tmp/rebuilt_500.png
magick compare -metric AE /tmp/live_500.png /tmp/rebuilt_500.png /tmp/diff_500.png; echo
```

## Azimuthal-current comet movie

`make_azimuthal_current_movie` reuses the three `2p_m*` runs and writes the
top-down 3-panel frames under `$TMPDIR/azim_movie_frames/azim_3panel_2d`.
Assemble them into a GIF with
[ImageMagick](https://imagemagick.org) (`brew install imagemagick`):

```bash
matlab -batch 'addpath("functions"); make_azimuthal_current_movie'
magick -delay 6 -loop 0 "$TMPDIR/azim_movie_frames/azim_3panel_2d/frame_"*.png movies/2p_pm1_azimuthal_current_3panel.gif
```

# Generating the manuscript figures

All figures are written to `figures/`. Copy them next to the manuscript `.tex`
before compiling. The table maps each figure to the runs it needs and the driver
that builds it:

| Fig. | Content | Required run(s) | Driver |
|------|---------|-----------------|--------|
| 1 | 1s electron on a sphere | `1s0_1fs` | `h_atom_post_processing` |
| 2 | 3D electron cloud (1s, 2p0) | `1s0_1ps`, `2p0_10ps` | `h_atom_post_processing` (2p0 cloud is full-resolution, ~12 min) |
| 3 | Distributions, 1s | `1s0_1ps` | `h_atom_post_processing` |
| 4 | Distributions, 2p0 | `2p0_10ps` | `h_atom_post_processing` |
| 5 | Distributions, (2,1,1) | `2p_m1_1ps` | `h_atom_post_processing` |
| 6 | Polar convergence, 2p0 | `2p0_10ps` | `h_atom_post_processing` |
| 7 | Radial convergence, 2s0 | `2s0_10ps` | `h_atom_post_processing` |
| 8 | Distribution deviations | `2p0_10ps`, `2s0_10ps` | `h_atom_post_processing` |
| 9 | Azimuthal circulation | `2p_m1_1ps`, `2p_m0_1ps`, `2p_mn1_1ps` | `make_azimuthal_current_figure` |
| 10 | Energies vs time | `2p0_10ps`, `2s0_10ps` | `h_atom_post_processing` |
| 11 | Kinetic-energy cutoff scan | `2s0`/`2p0` cutoff sweep (64 runs) | `analyse_cutoff_scan` |

Once the required runs are present under `data/`:

```bash
# Figs 1-8, 10 (needs 1s0_1fs, 1s0_1ps, 2p0_10ps, 2p_m1_1ps, 2s0_10ps).
# Includes Fig 2b (the full-resolution 2p0 cloud), which alone adds ~12 min.
matlab -batch 'addpath("functions"); h_atom_post_processing'

# Fig 9 (needs 2p_m1_1ps, 2p_m0_1ps, 2p_mn1_1ps)
matlab -batch 'addpath("functions"); make_azimuthal_current_figure'

# Fig 11 -- both panels: 2s0 = (2,0,0) left, 2p0 = (2,1,0) right (needs the cutoff sweeps above)
matlab -batch "addpath('functions'); analyse_cutoff_scan('2s0','figures','data','default');" 
matlab -batch "addpath('functions'); analyse_cutoff_scan('2p0','figures','data','default');"
```

By default `h_atom_post_processing` writes only the figures listed above. Set
`extended = true` to also generate diagnostic `(2,1,m)` figures that are not used
in the manuscript (per-state energies and xy clouds for the `m = 0, +1, -1`
states):

```bash
matlab -batch 'addpath("functions"); extended=true; h_atom_post_processing'
```

`analyse_cutoff_scan` pools every scan run found under `data/` (one curve per
`dt`, fixed y-range `[2, 9]x1e-19` J, analytic value as a dash-dot line) and
writes `figures/<state>_cutoff_scan_kinetic_manuscript.pdf`. Pass `'extended'`
as the 4th argument to get a 6-panel diagnostic figure instead.

```bash
matlab -batch "addpath('functions'); analyse_cutoff_scan('2s0','figures','data','extended')"
```
