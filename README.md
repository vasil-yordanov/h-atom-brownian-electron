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
params_set_name='2p1_1ps';
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

## $1s_0$ state
[![1s0_state](movies/1s0_1ps.gif)](movies/1s0_1ps.gif)

## $2p_0$ state
[![2p0_state](movies/2p0_10ps.gif)](movies/2p0_10ps.gif)

## $2p_1$ state
[![2p1_state](movies/2p1_1ps.gif)](movies/2p1_1ps.gif)


## $2s_0$ state
[![2s0_state](movies/2s0_10ps.gif)](movies/2s0_10ps.gif)
