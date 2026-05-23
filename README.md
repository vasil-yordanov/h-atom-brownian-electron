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

# Post processing
The simulation will creata 'data' directory for every parameters set. 
The post processing tooks around 5 min.
```Matlab
h_atom_post_processing
```

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
