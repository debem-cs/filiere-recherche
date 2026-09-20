# Nonlinear Kapitza pendulum simulation study

This folder separates numerical result generation from plotting and analysis:

- `KP_generate_results.m`: runs `KP_sim.slx` and writes the numerical results.
- `kP_grid_analysis.m`: loads existing results, calculates derived quantities
  and creates the static figures without running Simulink.
- `KP_DFT_view.m`: interactively displays a stored DFT after the user clicks
  a stable point in the amplitude-frequency grid.
- `KP_sim.slx`: contains only the nonlinear physical model and signal logging.

The purpose of this study is to complement the linear Floquet stability map
with nonlinear time-domain simulations. The stable pairs are not rediscovered
by simulation: they are read from `generated/kapitza_stability_data.mat`.

## 1. Physical model

The vertical displacement of the pendulum support is

```text
y(t) = A cos(omega t)
```

where:

- `A` is the displacement amplitude in metres;
- `omega = 2*pi*f` is the excitation angular frequency in rad/s;
- `f` is the excitation frequency in hertz.

The acceleration amplitude is therefore

```text
a0 = A*omega^2
```

The model uses the upright position as `theta = 0`. Its nonlinear equation is

```text
thetaDDot + 2*zeta*Omega*thetaDot
          + (-Omega^2 + a0*cos(omega*t)/L)*sin(theta) = 0
```

with

```text
Omega = sqrt(g/L).
```

The state vector is

```text
x = [theta; thetaDot].
```

Angles inside the calculation are expressed in radians.

## 2. Required Simulink interface

The MATLAB Function block should have four inputs:

```text
x, t, A, frequencyHz
```

Two Constant blocks can supply `A` and `frequencyHz`. Their values should be
the variable names `A` and `frequencyHz`, rather than fixed numbers. The study
script changes these variables before each simulation.

The Integrator initial-condition field must contain:

```text
x0
```

The To Workspace block must be connected to the complete state vector and use:

```text
Variable name: x
Save format: Timeseries
```

The complete MATLAB Function block code is:

```matlab
function xDot = fcn(x,t,A,frequencyHz)

g = 9.81;       % gravitational acceleration [m/s^2]
L = 1.0;        % pendulum length [m]
zeta = 0.01;    % damping ratio [-]

Omega = sqrt(g/L);
omega = 2*pi*frequencyHz;
a0 = A*omega^2;

phi = cos(omega*t);

theta = x(1);
thetaDot = x(2);

thetaDDot = -2*zeta*Omega*thetaDot ...
    - (-Omega^2 + (a0/L)*phi)*sin(theta);

xDot = [thetaDot; thetaDDot];

end
```

## 3. Amplitude-frequency grid

The engineer selects the vectors `amplitudeValues` and `frequencyValues` at
the beginning of `KP_generate_results.m`.

The supplied grid covers the selected physical amplitude and frequency
intervals:

```matlab
amplitudeValues = linspace(0.25,0.70,23);
frequencyValues = linspace(1.00,6.00,61);
```

For each `(A,f)` pair, the script converts the physical parameters to

```text
alpha = -Omega^2/(2*pi*f)^2
beta  = A/L.
```

It interpolates `log10SpectralRadius` from the saved Floquet map. Only points
with `log10SpectralRadius < 0` are sent to Simulink. At each stable pair, all
the following initial positions are tested:

```text
theta(0) = 5, 15, 25, 35, 45, 55, 65, 75 and 85 degrees
thetaDot(0) = 0 rad/s.
```

Only positive angles are considered because the ideal equation is symmetric
about `theta = 0`. This symmetry would no longer be sufficient if asymmetric
forces, stops or nonlinear damping were added.

## 4. Simulation duration

Every case is simulated for 1200 excitation periods:

```text
simulation time = 1200/f seconds.
```

Using a number of periods instead of one fixed physical time gives every
frequency the same number of forcing cycles.

The maximum solver step is one fiftieth of the excitation period. The model
may still use a variable-step solver. The final signal is resampled uniformly
before calculating its DFT.

## 5. Convergence and settling time

The upright equilibrium is

```text
theta = 0 modulo 2*pi
thetaDot = 0.
```

The angle is wrapped to the interval `[-pi,pi]`. The following equivalent
angular amplitude is then calculated:

```text
r = sqrt(theta^2 + (thetaDot/omega)^2).
```

Dividing velocity by the excitation frequency converts it to an equivalent
angular amplitude. This prevents the fast periodic velocity from making the
settling criterion unnecessarily restrictive. A trajectory is inside the
selected neighbourhood when `r <= 0.5 degree`.

It is classified as converged only if it remains inside that neighbourhood
for at least ten complete excitation periods before the simulation finishes.

Settling times are stored for every tested initial angle in
`settlingTimeAllAngles`. The displayed settling-time map uses the intermediate
initial condition of 35 degrees. This emphasizes nonlinear transient behaviour
instead of reproducing mainly the small-angle damping rate.

If 35 degrees is outside the attraction region at a Floquet-stable `(A,f)`
point, the displayed settling time is `NaN`. This means that the local upright
equilibrium is stable but that this particular initial condition did not reach
it during the simulation.

The settling time is not generally identical for all initial positions. The
local decay rate close to the equilibrium can be similar, but the preceding
nonlinear transient depends on the initial position. Large initial positions
can also converge to a different attractor or fail to converge.

A `NaN` value in `settlingTime` means that the 35-degree trajectory did not
satisfy the convergence criterion during the available simulation time. These
points are shown in grey, rather than being displayed as zero settling time.

## 6. Attraction result

The logical array `converged` stores the result for every amplitude, frequency
and tested initial angle.

The two-dimensional result `attractionAngle` contains the largest consecutive
tested angle that converged, starting at five degrees. For example:

```text
attractionAngle = 45 degrees
```

means that the simulations starting at 5, 15, 25, 35 and 45 degrees converged.
It does not claim that every continuous initial condition below 45 degrees
was tested.

This is therefore a practical angular estimate of the attraction region for
zero initial angular velocity. A complete region of attraction would also
require a sweep over initial angular velocity.

## 7. Discrete Fourier transform

The DFT remains calculated only for the five-degree simulation at every `(A,f)`
point inside the Floquet-stable region. It uses the final 32 excitation
periods. Points outside the Floquet-stable region contain `NaN` and are not
simulated.

The requested amplitude-frequency interval is preserved for the calculation.
For presentation, the result figures are automatically cropped to the
smallest rectangular interval containing all Floquet-stable grid points. This
removes unused grey margins without hiding unstable holes or nonconvergent
initial conditions located inside the stable interval.

The samples are first placed on a uniform time grid. The mean angle is then
removed, and MATLAB's `fft` is applied. An integer number of excitation periods
is used to limit spectral leakage and to align the forcing frequency with a
DFT bin.

All generated files are written to the `generated` folder beside the scripts.
The spectral results are saved separately in `generated/KP_DFT_results.mat`.
The saved quantities are:

- `dftFrequencyHz(:,frequencyIndex)`: physical frequency axis in hertz;
- `dftComplex(:,amplitudeIndex,frequencyIndex)`: complex positive-frequency
  coefficients, retaining spectral phase;
- `dftAmplitude(:,amplitudeIndex,frequencyIndex)`: one-sided angle amplitude
  spectrum in radians.

For a given grid point, the phase spectrum is obtained with:

```matlab
phaseRadians = angle(dftComplex(:,amplitudeIndex,frequencyIndex));
```

To inspect a result, run `KP_DFT_view.m` and click a green point in its
amplitude-frequency grid. The amplitude and phase plots are updated after
each click. Click another point to compare it, or press Enter in the grid
figure to finish.

The DFT window begins and ends at the same phase of the imposed cosine. Thus,
the phase at the forcing-frequency bin can be interpreted relative to the
base excitation. Phase should not be interpreted when the corresponding
spectral amplitude is approximately zero.

If the trajectory truly converges to the static upright equilibrium, its final
spectrum tends to zero. A nonzero spectrum indicates a remaining periodic,
subharmonic or non-periodic response rather than an exact static equilibrium.

## 8. Saved outputs

`KP_generate_results.m` creates exactly the two numerical files used by the
analysis:

```text
generated/KP_grid_results.mat
generated/KP_DFT_results.mat
```

`kP_grid_analysis.m` reads those files and creates:

```text
generated/KP_analysis_results.mat
generated/KP_settling_time.png
generated/KP_attraction_region.png
generated/KP_system_oscillation_energy.png
generated/KP_control_signal_energy.png
generated/KP_relevant_frequency_count.png
generated/KP_relevant_frequency_ratios.png
generated/KP_secondary_peak_strength.png
```

The generation files preserve the grid values, convergence classifications,
settling times, attraction-angle estimate and complete DFT data. The analysis
file contains only quantities derived from those existing results.

## 9. Energy maps

The positive residual oscillation-energy measure is calculated from the final
angle DFT using Parseval's relation:

```text
e_osc = 0.5*L^2*mean(thetaDot^2) + 0.5*g*L*mean(theta^2).
```

Because pendulum mass is not specified, this is specific energy in J/kg. It
describes residual oscillation in the final DFT window, not the signed
mechanical energy about the inverted equilibrium.

The control input is acceleration:

```text
u(t) = a0*cos(omega*t),  a0 = A*omega^2.
```

The plotted control quantity is its signal energy during one period:

```text
J_u = integral(u^2 dt) = 0.5*a0^2*T.
```

It is an acceleration-effort measure with units `(m/s^2)^2*s`, not actuator
energy in joules. Both energy maps use base-10 logarithmic colour values.

## 10. Relevant-frequency analysis

A frequency is counted when it is a local DFT maximum. The strongest peak
`a_1` is always retained. Every additional peak must satisfy:

```text
a_k > 0.10*a_1
```

There is no absolute amplitude threshold. Consequently, every Floquet-stable
point with stored DFT data has at least one relevant frequency. A point may
still have extremely small physical oscillation amplitude; that information
is shown separately by the residual oscillation-energy map.

Here `a_1` means the largest peak of the pendulum-output spectrum. It is not
the amplitude of the control input and its frequency `f_1` does not have to be
equal to the excitation frequency `f`.

Only frequencies between zero and five times the excitation frequency are
examined. Candidate peaks separated by no more than two DFT bins are treated
as one spectral component. This avoids counting adjacent leakage bins as
different physical frequencies.

The analysis stores the number of relevant frequencies and the three strongest
accepted peaks. Frequency maps use the normalized ratios `f_1/f`, `f_2/f` and
`f_3/f`. Separate maps show `a_2/a_1` and `a_3/a_1`, so a second or third peak
that is nearly as strong as the first remains visible.

The three frequency-ratio panels use the same colour limits. Their upper limit
is the largest stored frequency ratio rounded upward to the next integer; it
is not fixed at the five-harmonic search limit.

## 11. Interpreting the final DFT

The first useful spectral results are two scalar maps:

1. Response amplitude at the excitation frequency `f`.
2. Response phase at `f`, displayed only where that amplitude is significant.

The next quantities to inspect are the amplitudes at `f/2` and `2*f`. A strong
component at `f/2` indicates a period-two response. A broad spectrum can point
to a non-periodic response, but it should be confirmed from the time history
and a longer simulation.
