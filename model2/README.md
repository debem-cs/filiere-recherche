# Model 2: viscously damped pendulum and feedback control

This folder contains the second nonlinear pendulum model and the closed-loop
models used in the project. Unlike Model 1, Model 2 represents dissipation as
viscous angular damping and does not include a base-velocity damping term.

The folder supports two activities:

1. open-loop simulation under a prescribed vertical excitation;
2. comparison of the implemented vibrational feedback controller with a
   non-vibrational `LgV` controller.

## Open-loop physical model

Define

```math
\Omega=\sqrt{\frac{g}{L}},\qquad
u(t)=A\omega^2\cos(\omega t),\qquad
\omega=2\pi f.
```

The nonlinear equation is

```math
\ddot\theta+2\zeta\Omega\dot\theta+
\left[-\Omega^2+\frac{u(t)}{L}\right]\sin\theta=0.
```

Equivalently,

```math
\ddot\theta=-2\zeta\Omega\dot\theta+
\left[\frac{g}{L}-\frac{A\omega^2}{L}\cos(\omega t)\right]\sin\theta.
```

The upright equilibrium is `theta = 0`. The open-loop parameter vector is

```matlab
p = [g; L; zeta; A; f];
```

and the initial state is

```matlab
x0 = [theta0; thetaDot0];
```

## Files

### `kapitza_manual_model2.slx`

Open-loop Simulink model. It logs:

- `theta_ts`: angle in degrees;
- `omega_ts`: angular velocity in rad/s.

### `parameters_kapitza_model2.m`

Configures and runs the open-loop model. The currently saved case uses

```text
g      = 9.81 m/s^2
L      = 0.25 m
lambda = 1 s^-1
zeta   = lambda/(2*sqrt(g/L))
A      = 0.08 m
f      = 15 Hz
x0     = [45 deg; 0 rad/s]
```

The script uses `ode45`, relative tolerance `1e-7`, absolute tolerance `1e-9`,
maximum step `1/(100*f)`, and a 20 s simulation horizon.

### `model2_control.slx`

Closed-loop model with the implemented vibrational input

```math
u_{\mathrm{vib}}(t)=
A\omega^2\cos(\omega t)-\sqrt{\omega}V_2\sin(\omega t),
```

where

```math
V_2=\frac{1}{2}(c_2\dot\theta+\theta)^2,
\qquad
A=\frac{2\alpha_c}{\omega^{3/2}}.
```

The implementation uses `theta` inside `V2`, rather than `sin(theta)`. This is
the modified controller evaluated in the project and should not be interpreted
as an unchanged implementation of the reference control law.

### `model2_control_LgV.slx`

Comparison model with the non-vibrational input

```math
u_{\mathrm{LgV}}=
\frac{\alpha_c c_2}{L}\sin\theta
\left(c_2\dot\theta+\sin\theta\right).
```

### `kapitza_controle.m`

Runs `model2_control.slx` and `model2_control_LgV.slx` with identical plant,
initial-condition, solver, and horizon settings, then plots the two responses.
The script defines

```matlab
p  = [g; L; zeta; A; f];
pc = [c1; c2; alpha_ctrl];
```

with the saved controller parameters

```text
omega      = 500 rad/s
alpha_ctrl = 40
c1         = 1
c2         = 1.5
x0         = [150 deg; 0 rad/s]
```

Both models use `ode45`, relative tolerance `1e-7`, absolute tolerance `1e-9`,
maximum step `1/(100*f)`, and a 10 s simulation horizon. The angle and angular
velocity are available through `theta_ts` and `omega_ts`; the vibrational model
also provides the logged signal `Vx` used by the script.

## Running the simulations

Open MATLAB in this directory, or make this directory the current MATLAB
folder.

For the open-loop simulation:

```matlab
parameters_kapitza_model2
```

For the controller comparison:

```matlab
kapitza_controle
```

## Changing a test case

- Change `A`, `f`, or `x0` in `parameters_kapitza_model2.m` for an open-loop
  experiment.
- Change `w`, `alpha_ctrl`, `c1`, `c2`, or `x0` in `kapitza_controle.m` for a
  feedback experiment.
- Keep angles in radians in `x0`, even though the logged angle and plot labels
  are expressed in degrees.
- Preserve a maximum solver step that resolves the high-frequency carrier.
