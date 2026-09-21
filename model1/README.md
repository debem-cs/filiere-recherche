# Model 1: pendulum with base-velocity damping contribution

This folder contains the first nonlinear model used in the project. It follows
a damping representation in which the force depends on the relative motion of
the pendulum mass and the vertically moving base.

The model is intended for open-loop studies of how excitation amplitude,
frequency, and initial angle affect stabilization of the upright equilibrium.

## Equation of motion

The support displacement, velocity, and acceleration are prescribed as

```math
z_b(t)=-A\cos(\omega t),\qquad
v_b(t)=A\omega\sin(\omega t),\qquad
a_b(t)=A\omega^2\cos(\omega t),
```

with `omega = 2*pi*f`. The implemented nonlinear equation is

```math
\ddot\theta=
\frac{g}{L}\sin\theta
-\frac{b}{m}\dot\theta
+\frac{b}{mL}v_b(t)\sin\theta
-\frac{a_b(t)}{L}\sin\theta.
```

Here:

- `theta = 0` is the upright equilibrium;
- `L` is the pendulum length;
- `m` is the point mass;
- `b` is a dimensional damping coefficient;
- `A` is the base-displacement amplitude;
- `f` is the excitation frequency in hertz.

The parameter `b` is not a dimensionless damping ratio. This model differs
from Model 2 because it includes the term proportional to the base velocity
`v_b(t)`.

## Files

### `kapitza_manual.slx`

Simulink implementation of the nonlinear equation. The model receives the
state, time, and parameter vector, integrates the state derivative, and logs:

- `theta_ts`: angle in degrees;
- `omega_ts`: angular velocity in rad/s.

The Integrator initial condition is the workspace variable `x0`.

### `parameters_kapitza.m`

Defines one simulation case, loads `kapitza_manual.slx`, runs the model, and
plots the angle and angular velocity. The parameter vector passed to Simulink
is

```matlab
p = [g; l; m; b; A; f];
```

and the initial state is

```matlab
x0 = [theta0; thetaDot0];
```

where both state entries use SI angular units: radians and rad/s.

The values currently saved in the script are:

```text
g  = 9.81 m/s^2
l  = 0.25 m
m  = 1.0 kg
b  = 0.0628 kg/s
A  = 0.04 m
f  = 15 Hz
x0 = [75 deg; 0 rad/s]
```

### `exp_kapitza.m`

Runs four amplitude--frequency cases and displays their angle trajectories in
a tiled figure:

```text
A = 0.04 m, f =  7 Hz
A = 0.04 m, f = 15 Hz
A = 0.02 m, f = 15 Hz
A = 0.04 m, f = 20 Hz
```

For each case, the script sets a 100 s simulation horizon and limits the
maximum step to `1/(100*f)`. The complete numerical outputs are retained in the
workspace cell array `resultados`.

## Running the simulations

Open MATLAB in this directory, or make this directory the current MATLAB
folder.

For the single default case:

```matlab
parameters_kapitza
```

For the four-case comparison:

```matlab
exp_kapitza
```

`exp_kapitza.m` begins by running `parameters_kapitza.m`, so it first creates
the default-case figure and then executes the comparison sweep.

## Changing a test case

Edit `A`, `f`, or `x0` in `parameters_kapitza.m` for a single run. For the
comparison study, edit the rows of `casos` in `exp_kapitza.m`. When increasing
the excitation frequency, keep the maximum solver step small relative to the
excitation period.
