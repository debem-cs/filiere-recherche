# Robust stabilization of a parametrically excited pendulum

This repository contains the MATLAB and Simulink material developed for the
Group 30 research project on robust stabilization of parametrically excited
dynamical systems. The Stephenson--Kapitza pendulum is used as the main
benchmark.

The project covers:

- two nonlinear models of a vertically excited inverted pendulum;
- direct Floquet stability-map construction;
- nonlinear simulations selected from the Floquet-stable region;
- time-domain and frequency-domain analysis of the nonlinear response;
- comparison of vibrational and non-vibrational feedback controllers.

The angular convention is the same throughout the repository: `theta = 0`
denotes the upright equilibrium, and positive angles are measured from the
upper vertical. Internal angular calculations use radians unless stated
otherwise.

## Repository structure

```text
MATLAB/
|-- model1/
|   |-- kapitza_manual.slx
|   |-- parameters_kapitza.m
|   |-- exp_kapitza.m
|   `-- README.md
|-- model2/
|   |-- kapitza_manual_model2.slx
|   |-- model2_control.slx
|   |-- model2_control_LgV.slx
|   |-- parameters_kapitza_model2.m
|   |-- kapitza_controle.m
|   `-- README.md
|-- floquet-stability/
|   |-- stability-map/
|   `-- sim-grid-analysis/
`-- README.md
```

### `model1`

Implements the nonlinear pendulum model with a damping force that contains
both angular-velocity and base-velocity contributions. It is used for the
initial amplitude--frequency experiments.

See [`model1/README.md`](model1/README.md).

### `model2`

Implements the nonlinear model with viscous angular damping. The folder also
contains the closed-loop Simulink models used to compare the implemented
vibrational controller with the non-vibrational `LgV` controller.

See [`model2/README.md`](model2/README.md).

### `floquet-stability/stability-map`

Constructs the local stability map of the upright equilibrium by integrating
the periodic linearized system over one excitation period and evaluating the
spectral radius of the monodromy matrix.

See
[`floquet-stability/stability-map/README.md`](floquet-stability/stability-map/README.md).

### `floquet-stability/sim-grid-analysis`

Uses the stable points from the Floquet map as inputs to nonlinear Simulink
experiments. It evaluates finite-angle recovery, settling time, residual
oscillation, acceleration-signal effort, and frequency content.

See
[`floquet-stability/sim-grid-analysis/README.md`](floquet-stability/sim-grid-analysis/README.md).

## Requirements

- MATLAB;
- Simulink;
- Symbolic Math Toolbox for the symbolic linearization used by the
  Floquet-map script.

No installation step is required. Clone the repository and run each script
from the directory in which it is stored so that MATLAB can find the associated
Simulink models and generated data.

## Quick start

### Model 1: open-loop simulations

```matlab
cd model1
parameters_kapitza
```

To reproduce the four amplitude--frequency comparisons:

```matlab
exp_kapitza
```

### Model 2: open-loop simulation

```matlab
cd model2
parameters_kapitza_model2
```

### Model 2: feedback comparison

```matlab
cd model2
kapitza_controle
```

### Floquet stability map

```matlab
cd(fullfile('floquet-stability','stability-map'))
kapitza_stability_map
```

The calculation writes its data and figures to the local `generated/`
directory.

### Nonlinear grid analysis

The nonlinear grid study has its own execution order and data requirements.
Follow the instructions in
[`floquet-stability/sim-grid-analysis/README.md`](floquet-stability/sim-grid-analysis/README.md)
before running its scripts.

## Numerical-resolution note

The excitation can be much faster than the pendulum response. The scripts that
configure the solver limit the maximum integration step to a fraction of the
excitation period. If the parameters or solver configuration are changed,
retain a sufficiently small maximum step; otherwise, skipped excitation cycles
can produce misleading apparent stability.

## Repository

The project is hosted at
<https://github.com/debem-cs/filiere-recherche>.
