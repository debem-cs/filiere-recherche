# Damped Kapitza pendulum: explicit Floquet stability map

This project calculates the local Floquet stability of the upward equilibrium
of a vertically excited, damped Kapitza pendulum. It deliberately uses a direct
engineering-style implementation: the physical parameters, differential
equation, parameter grid, time discretization, and stability calculation appear
in their mathematical order.

The project contains two MATLAB files:

```text
kapitza_floquet_explicit/
|-- kapitza_stability_map.m   equations and numerical calculation
|-- plot_stability_maps.m     plotting and figure formatting only
|-- generated/                results created when the script runs
`-- README.md
```

Run `kapitza_stability_map.m`. Symbolic Math Toolbox is required for the
symbolic definition and linearization of the nonlinear equation.

## 1. Physical model

The pendulum has length `L`, gravitational acceleration `g`, and physical
damping ratio `zeta`. Define

```math
\Omega=\sqrt{\frac{g}{L}}.
```

The angle `theta = 0` represents the upward equilibrium. The vertical
acceleration of the suspension is prescribed as

```math
u(t)=a_0\phi(\omega t),
```

where:

- `a0` is acceleration amplitude in `m/s^2`;
- `omega` is excitation angular frequency in `rad/s`;
- `phi` is a dimensionless periodic waveform;
- `phi(tau + 2*pi) = phi(tau)` in the present implementation.

The nonlinear damped equation is

```math
\ddot\theta+2\zeta\Omega\dot\theta+
\left[-\Omega^2+\frac{a_0}{L}\phi(\omega t)\right]\sin\theta=0.
```

The sign of the excitation depends on the chosen vertical coordinate and input
convention. For a sinusoid, changing the sign is equivalent to shifting the
input phase by half a period and therefore does not change the Floquet
multipliers.

## 2. Relation to base-displacement amplitude

The physical frequency-amplitude figure uses `A` for base-displacement
amplitude. For sinusoidal support motion, acceleration amplitude and
displacement amplitude satisfy

```math
a_0=A\omega^2.
```

For example, the support displacement

```math
y_b(t)=-A\cos(\omega t)
```

has acceleration

```math
\ddot y_b(t)=A\omega^2\cos(\omega t).
```

Thus, with `phi(tau) = cos(tau)`, the acceleration input used in the pendulum
equation is

```math
u(t)=A\omega^2\cos(\omega t).
```

The symbols `a0` and `A` are not two independent inputs: `a0` is acceleration
amplitude and `A = a0/omega^2` is the corresponding displacement amplitude.

## 3. Nondimensional time and parameters

Introduce the phase variable

```math
\tau=\omega t.
```

The derivative transformations are

```math
\frac{d}{dt}=\omega\frac{d}{d\tau},\qquad
\frac{d^2}{dt^2}=\omega^2\frac{d^2}{d\tau^2}.
```

Dividing the physical equation by `omega^2` gives

```math
\theta''+2\delta\theta'+
\left[\alpha+\beta\phi(\tau)\right]\sin\theta=0,
```

where primes denote derivatives with respect to `tau`, and

```math
\alpha=-\frac{\Omega^2}{\omega^2},\qquad
\beta=\frac{a_0}{L\omega^2},\qquad
\delta=\frac{\zeta\Omega}{\omega}.
```

For sinusoidal base displacement, `a0 = A*omega^2`, hence

```math
\beta=\frac{A}{L}.
```

This is why a rectangular `(alpha,beta)` diagram maps naturally to frequency
versus base-displacement amplitude.

For the physical upright-pendulum branch, `alpha < 0`, and

```math
\omega=\frac{\Omega}{\sqrt{-\alpha}},\qquad
\delta=\zeta\sqrt{-\alpha}.
```

The code also permits `alpha > 0` to display the conventional stable-oscillator
extension of the Hill diagram. It uses `delta = zeta*sqrt(abs(alpha))` on that
extension. Positive `alpha` does not map to the physical upright Kapitza
pendulum and is therefore excluded from the physical frequency plot.

## 4. Symbolic nonlinear model and linearization

Using the normalized state

```math
x=\begin{bmatrix}\theta & \theta'\end{bmatrix}^{T},
```

the nonlinear first-order system is

```math
x'=
\begin{bmatrix}
\theta'\\
-2\delta\theta'-(\alpha+\beta\phi(\tau))\sin\theta
\end{bmatrix}.
```

The script defines this equation symbolically and computes its Jacobian:

```math
A(\tau,x)=\frac{\partial f}{\partial x}.
```

At the upward equilibrium `theta = theta' = 0`, the periodic linearization is

```math
x'=A(\tau)x,
\qquad
A(\tau)=
\begin{bmatrix}
0 & 1\\
-\alpha-\beta\phi(\tau) & -2\delta
\end{bmatrix}.
```

Equivalently, the scalar linear equation is

```math
\theta''+2\delta\theta'+
(\alpha+\beta\phi(\tau))\theta=0.
```

Floquet stability of this equation determines the local stability of the
upward equilibrium of the nonlinear system.

## 5. State-transition and monodromy matrices

For a linear time-periodic system

```math
x'=A(\tau)x,
```

the state-transition matrix satisfies

```math
\frac{d\Phi}{d\tau}=A(\tau)\Phi,
\qquad
\Phi(0)=I.
```

It maps an initial perturbation to a later perturbation:

```math
x(\tau)=\Phi(\tau,0)x(0).
```

The monodromy matrix is the transition over one complete forcing period:

```math
M=\Phi(2\pi,0).
```

After `k` complete periods,

```math
x(2\pi k)=M^k x(0).
```

The eigenvalues of `M` are the Floquet multipliers.

## 6. Explicit time discretization

The interval `[0,2*pi]` is divided into `N` equal intervals:

```math
h=\frac{2\pi}{N},\qquad
\tau_k=kh.
```

The coefficient matrix is sampled at each midpoint:

```math
\tau_{k+1/2}=\left(k+\frac12\right)h.
```

During that short interval, `A(tau)` is approximated as constant. The local
transition is therefore

```math
E_k=\exp\left[A(\tau_{k+1/2})h\right].
```

The complete transition is assembled chronologically:

```math
M=E_{N-1}E_{N-2}\cdots E_1E_0.
```

In the code, this is implemented by left multiplication. Reversing the order
would generally give an incorrect result because coefficient matrices at
different times do not necessarily commute.

Midpoint freezing is second-order accurate for a smooth coefficient matrix.
It is more accurate than freezing the matrix at the beginning of each interval,
which is only first-order accurate.

## 7. Analytical exponential of each frozen matrix

Calling MATLAB's general `expm` function millions of times is unnecessarily
expensive for this `2-by-2` system. At one midpoint, write

```math
q_k=\alpha+\beta\phi(\tau_{k+1/2}),
```

so that

```math
A_k=
\begin{bmatrix}
0 & 1\\
-q_k & -2\delta
\end{bmatrix}
=-\delta I+B_k.
```

The remaining matrix satisfies

```math
B_k^2=(\delta^2-q_k)I.
```

Define

```math
D_k=\delta^2-q_k.
```

If `Dk > 0`, the local transition is evaluated using hyperbolic functions:

```math
E_k=e^{-\delta h}
\left[
\cosh(\sqrt{D_k}h)I+
\frac{\sinh(\sqrt{D_k}h)}{\sqrt{D_k}}B_k
\right].
```

If `Dk < 0`, let `nu = sqrt(-Dk)`. The equivalent real expression is

```math
E_k=e^{-\delta h}
\left[
\cos(\nu h)I+
\frac{\sin(\nu h)}{\nu}B_k
\right].
```

At `Dk = 0`, the limiting expression is

```math
E_k=e^{-\delta h}(I+hB_k).
```

These expressions are the exact matrix exponential for the frozen interval;
the approximation arises only from treating `A(tau)` as constant during that
interval.

## 8. Floquet stability criterion

Let the Floquet multipliers be

```math
\mu_i=\operatorname{eig}(M).
```

Define the spectral radius

```math
\rho(M)=\max_i|\mu_i|.
```

For the damped externally forced system:

- `rho(M) < 1`: locally exponentially stable;
- `rho(M) > 1`: unstable;
- `rho(M) = 1`: transition boundary.

The plotted quantity is

```math
S(\alpha,\beta)=\log_{10}\rho(M).
```

Therefore:

- `S < 0` is stable;
- `S > 0` is unstable;
- the black contour `S = 0` is the Floquet stability boundary.

## 9. Determinant verification

Liouville's formula provides an independent numerical check:

```math
\det M=\exp\left(\int_0^{2\pi}\operatorname{tr}A(\tau)\,d\tau\right).
```

Because

```math
\operatorname{tr}A(\tau)=-2\delta,
```

the exact determinant is

```math
\det M=e^{-4\pi\delta}.
```

The script compares the determinant of every computed monodromy matrix with
this result and reports the maximum relative error.

This check can reveal errors in transition ordering, time scaling, or the local
matrix-exponential formula. A small determinant error does not by itself prove
that the time grid is sufficiently fine, so convergence must also be checked.

## 10. Parameter grids and convergence

The current requested map is

```matlab
alphaValues = linspace(-1.0,0.25,1000);
betaValues = linspace(0.0,4.0,300);
numberOfTimeSteps = 120;
```

There are two independent discretizations:

1. The `(alpha,beta)` grid controls the spatial resolution of the plotted map.
2. `numberOfTimeSteps` controls the accuracy of each monodromy matrix.

Increasing the number of alpha and beta points makes the plotted boundary
smoother but does not improve the time propagation at an individual point.
Increasing `numberOfTimeSteps` improves the monodromy approximation but does
not increase the visual grid density.

A practical convergence check is:

1. Calculate the map with `N` time steps.
2. Repeat with `2*N` time steps.
3. Compare `log10SpectralRadius` and the location of the zero contour.
4. Pay particular attention to points near `rho(M) = 1`, where a small
   numerical error can change the stable/unstable classification.

For a discontinuous waveform, switching phases should coincide exactly with
interval boundaries. A step should not cross a discontinuity.

## 11. Physical frequency-amplitude conversion

Only the `alpha < 0` portion corresponds to the upward Kapitza pendulum. Its
physical angular frequency and ordinary frequency are

```math
\omega=\frac{\Omega}{\sqrt{-\alpha}},\qquad
f=\frac{\omega}{2\pi}.
```

For sinusoidal base motion, displacement amplitude is

```math
A=\beta L.
```

The plotting function uses these formulas to create the dimensional frequency
versus displacement-amplitude figure. The amplitude unit is selected
automatically as metres, centimetres, millimetres, or micrometres.

`maximumFrequencyRatioShown` limits `omega/Omega` in the dimensional plot.
Setting it to `Inf` includes every negative-alpha grid point. High frequencies
correspond to alpha values very close to zero from below, so a nonuniform alpha
grid may be useful when high-frequency resolution is important.

## 12. Changing the waveform

Edit these two lines in `kapitza_stability_map.m`:

```matlab
phi = cos(tau);
phiName = 'cos(\tau)';
```

For example:

```matlab
phi = cos(tau) + sym(1)/5*cos(2*tau);
phiName = 'cos(\tau)+0.2cos(2\tau)';
```

The waveform samples are evaluated symbolically once and reused throughout the
parameter sweep.

The optimized local transition relies on the scalar second-order structure

```math
\theta''+2\delta\theta'+q(\tau)\theta=0.
```

Changing `phi` preserves that structure and requires no other code changes.
Changing to a fundamentally different system requires deriving its new state
matrix and adapting `calculateMonodromyMatrix` accordingly.

## 13. Generated results

When saving is enabled, the `generated` folder contains:

- `kapitza_stability_data.mat`;
- `kapitza_stability_map_alpha_beta.png`;
- `kapitza_stability_map_frequency_amplitude.png`.

The MAT-file contains the parameter vectors, stability metric, logical stability
map, determinant-error diagnostic, physical constants, time-step count, and
waveform name.
