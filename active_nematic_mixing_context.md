# Project Context: Passive-Scalar Mixing in 2D Active Nematic Turbulence

This document is a self-contained briefing for an LLM assisting with a computational
physics project. It captures the model, the numerical setup, the theoretical scaling
framework, and the key results obtained so far. Read it before answering questions
about the project.

---

## 1. The physical system

A 2D **active nematic** fluid (director field + velocity field) with a **passive
scalar** (concentration field) advected by the nematic flow. The scalar is a circular
"drop" that is stirred by the self-generated active flow and simultaneously diffuses.
The scalar is **one-way coupled**: the flow stirs the scalar, but the scalar does not
feed back on the nematic.

The flow is in the **Stokes (inertialess, Re = 0) regime**: the velocity field has no
dynamics of its own; it is slaved instantaneously to the nematic stress. The velocity
is nonetheless strongly time-dependent because its *source* (the director field theta)
is a live dynamical field. One never time-steps the velocity — only theta and the
scalar; the velocity is recovered each step by inverting a biharmonic equation.

---

## 2. Governing equations

Stream-function / vorticity / director formulation. Symbols: psi = stream function,
omega = vorticity, theta = director angle, c = scalar concentration.

```
(6A)  Lap(omega) = -Lap^2(psi) = s(r,t)

(6B)  s(r,t) = (1/2)(R/A) Lap^2(theta)
             + S[ (1/2)(d_xx - d_yy) sin(2 theta) - d_xy cos(2 theta) ]

(7)   d_t theta + (d_y psi)(d_x theta) - (d_x psi)(d_y theta)
             + (1/2) Lap(psi) = (1/A) Lap(theta)

      d_t c + u . grad(c) = D Lap(c),   with u = (d_y psi, -d_x psi)
```

Velocity from stream function: **u = (d_y psi, -d_x psi)**.

### Term-by-term meaning
- `Lap^2(psi) = s` is **2D biharmonic Stokes flow**: u is slaved to the source s
  (curl of the internal force density).
- Active stress term `S[...]` is the **engine**: orientation gradients generate flow.
- `(R/A) Lap^2(theta)` is **elastic backflow** (director relaxation pushes fluid).
- In (7): `u.grad(theta)` = advection of orientation; `(1/2)Lap(psi) = -omega/2` =
  **corotation** (director rotates with local vorticity); `(1/A)Lap(theta)` =
  **rotational diffusion / elastic relaxation**.
- **No flow-alignment (strain-coupling) term** is present — this is a corotation-only
  model. Consequence: +1/2 defect self-propulsion is driven purely by active backflow,
  so defect motility may differ quantitatively from full Beris-Edwards nematics. State
  this caveat in any defect study.

---

## 3. Physical meaning of parameters (CRITICAL — non-standard convention)

In these **non-dimensionalized** equations:

- **S = +/-1 only.** Time is measured in units of the active time, which normalizes the
  active stress coefficient to unity, leaving only its sign. **S = +1 = extensile**
  (pushers, bend-unstable); **S = -1 = contractile** (pullers, splay-unstable). S sets
  the *symmetry*, NOT the magnitude of activity.
- **A is the activity / Ericksen number** (NOT merely a rotational viscosity). Because
  the activity magnitude was absorbed into the rescaling, A is the single knob that
  controls how unstable/turbulent the system is. 1/A is the orientational diffusivity;
  with activity fixed at 1, large A = weak relaxation relative to activity = MORE
  unstable, finer-scaled, more turbulent. **Large A = "more active."** Turning up A
  makes the active vortices smaller and more numerous (not faster).
- **R is the elastic-backflow strength** (dimensionless). It always enters as
  **(1 + R/4)** and is purely **stabilizing**: larger R raises the threshold, enlarges
  the active length, reduces defect density. It never drives anything.
- **D / Diff is the scalar diffusivity.** IMPORTANT: in the code the effective scalar
  diffusivity is **D_mol = Diff/A** (Diff is divided by A). So A enters the scalar
  problem too, tying molecular diffusion to activity.

---

## 4. Key scaling results (derived from linear stability of the equations)

Linear growth rate about a uniform director:
```
sigma(k) = (S/2) cos(2 phi)  -  (1 + R/4) k^2 / A
```
The destabilizing term is **k-independent** (saturates at S/2); damping grows as k^2.
Long wavelengths are always unstable in an infinite domain; a finite box sets threshold.

**Active length** (fundamental scale, size of active vortices/defect spacing):
```
l_a = sqrt( 2 (1 + R/4) / A )     ->   l_a ~ A^(-1/2)
```

**Velocity scale** (S = 1):  `u ~ l_a ~ A^(-1/2)`   (verified: measured u_max ≈ 0.014 at
A = 10000; predicts u(A) ∝ A^(-1/2)).

**Strain (stirring) rate:**  `gamma_dot ~ u / l_a ~ S ~ O(1)`  — independent of A.
So the dynamical timescale of the turbulence is O(1) in these units; raising A makes
structures finer, not faster.

**Spontaneous-flow threshold (box size L):**
```
S_c ~ (2 + R/2) / A * (2 pi / L)^2     i.e. onset at activity number A ~ (2 pi)^2 ~ 39.5
```

### Governing dimensionless groups
- **Activity number:**  `Activity = (L/l_a)^2 = A L^2 / [2(1+R/4)]`
- **Elastic ratio:** R (enters as 1 + R/4)
- **Peclet number (drop scale):**  `Pe = gamma_dot * r0^2 / D_mol = r0^2 A / Diff`
- **Geometric:** r0 / l_a (drop radius vs active length)
- **Batchelor scale:**  `l_B = sqrt(D_mol / gamma_dot) = sqrt(Diff/A)` ;
  `l_B / l_a = sqrt( Diff / [2(1+R/4)] )`  — **independent of A**, depends only on Diff.

### Timescales
- Active (stirring) time: `tau_a = 1/gamma_dot ~ 1/S ~ O(1)` — the natural clock.
- Diffusive time of the drop: `tau_D = r0^2 / (2 D_mol) = r0^2 A / (2 Diff)`.
- Ratio of clocks is fixed by the activity number; time carries no new dimensionless
  group. **Non-dimensionalize dynamics by tau_a (i.e. by 1/S).**

---

## 5. The eddy-diffusivity prediction (central result of the project)

When the drop is much larger than the active vortices, the turbulence transports the
scalar like an effective diffusion (mixing-length argument):
```
D_turb ~ kappa * u * l_a ~ kappa * l_a^2 = kappa * 2(1+R/4)/A    (proportional to 1/A)
```
Both u (~A^-1/2) and l_a (~A^-1/2) contribute a half-power, combining to 1/A.

Total effective diffusivity:
```
D_eff = D_mol + D_turb = [ Diff + const ] / A ,   const = 2 kappa (1+R/4)
```
**Both channels scale as 1/A.** This is why curves of c_max(t) at different activities
**collapse when time is divided by A** (t/A). The collapse is the experimental
fingerprint of this scaling. Note: the t/A collapse is partly *built into* the code by
the Diff/A convention — at fixed Diff, the active-scale Peclet l_B/l_a is A-independent,
so an A-sweep changes the clock speed but not the mixing *character*. To change the
stirring-vs-diffusion competition you must sweep **Diff**, not A.

---

## 6. Regime map (exponential vs power-law decay of the scalar maximum)

Reference framework: **Villermaux, "Mixing Versus Stirring," Annu. Rev. Fluid Mech.
2019, 51:245-273** (lamellar theory of mixing). Key relations from that paper:
- 1D lamella max concentration: `theta(tau) = erf(1/(4 sqrt(tau))) -> tau^(-1/2)` in
  *warped* (Ranz) time tau = D * integral dt'/s(t')^2 — NOT real time.
- Generalizes to d dimensions: `theta(t) ~ (s0/sqrt(Dt))^d`. For an **isotropic 2D
  blob in pure diffusion, d = 2 gives theta ~ 1/t** (NOT 1/sqrt(t); the 1/sqrt(t) is
  the 1D-lamella result).
- For chaotic/exponential stretching the real-time decay is **exponential** ~e^(-gamma_dot t)
  after a mixing time t_s ~ (1/(2 gamma_dot)) ln(Pe).

**Resulting regime map for THIS project:**

| Regime | Condition | c_max(t) behavior | Collapses in |
|---|---|---|---|
| Pure diffusion (no flow) | control | plateau then **t^(-1)** (2D) | t/A |
| Effective diffusion | Diff >~ 1  (l_B >~ l_a, no lamellar window) | plateau then **t^(-1)**, then exp box-mode tail | t/A |
| Lamellar (Villermaux) | Diff << 1 (l_B << l_a) | plateau to t_s then **~exp(-gamma_dot t)** | raw t |

At Diff ~ 1, l_B/l_a ~ 0.7 -> **no lamellar window exists**; the system is in the
effective-diffusion regime. Decisive experiment to confirm the map: **sweep Diff** (not
A). Curves should collapse in t/A at high Diff (power-law -1) and shift to raw-t
collapse with an exponential mid-stage at low Diff (~0.01).

Additional refinements from Villermaux:
- c_max (the single grid maximum) is carried by the **least-stretched** material element,
  i.e. the weak tail of the FTLE / stretching-rate distribution -> intermittent, jumpy.
  Mean/variance measures (scalar dissipation chi = D<|grad c|^2>, scalar variance) are
  far smoother diagnostics.
- Late-time **overlap/aggregation** (lamellae merging in the periodic box) eventually
  drives p(c) toward a Gamma distribution; c_max then relaxes toward the mean.

---

## 7. Numerical setup and the IMEX fix (IMPORTANT for stability)

- **Pseudo-spectral** solver (Julia), FFT-based, with **2/3 dealiasing**. Grid N x N,
  domain side L. Standard run: **N = 512, L = 1** (so Delta_x = 1/N ~ 0.002).
  Earlier runs used N = 128.
- The flow is recovered each step by inverting the biharmonic for psi from theta; the
  velocity is never time-stepped.
- **Original scheme was fully explicit (forward Euler)** for theta and the scalar. The
  diffusion stability limit is `dt < 2A / (Diff * max(k2))`, which is very restrictive
  at large N because max(k2) grows as N^2. At N = 512 with dt = 0.02 the explicit scheme
  BLOWS UP (worst at the *lowest* Diff, since effective diffusivity Diff/A is largest
  there ... actually the stiff linear term is theta's (1+R/4)/A and the scalar's Diff/A).
- **FIX (adopted): semi-implicit / IMEX Euler.** Treat the stiff linear diffusion
  implicitly (in Fourier space, just divide), keep advection + active driving explicit.
  This is unconditionally stable for the diffusion; only the advective CFL remains, which
  is loose because velocities are small (u_max ~ 0.014-0.02). Implemented updates:
  ```julia
  # theta update (active driving k2*psi/2 explicit, diffusion implicit):
  theta_new_tr = (theta_tr .+ dt .* (cov_der_the_tr .+ dealias .* k2 .* psi_tr ./ 2)) ./ (1 .+ dt .* k2 ./ A)

  # scalar update (advection explicit, diffusion implicit):
  covc_tr  = dealias .* fft(covc_temp)          # covc_temp = dypsi.*dxconc .- dxpsi.*dyconc
  conc_tr  = (conc_tr .- dt .* covc_tr) ./ (1 .+ dt .* (Diff/A) .* k2)
  conc_out = real.(ifft(conc_tr))
  ```
  (Apply the SAME implicit update to the second scalar field conc2.)
  Naming note: this is **semi-implicit (IMEX)**, not fully implicit Euler.
  A small residual elastic damping (the -R/(2A) Ericksen piece, ~ (R/4)k^2/A * theta)
  stays explicit; fine for R ~ O(1).
- With IMEX, **dt = 0.02 is fine across A = 5000-20000** (advective CFL allows ~0.05;
  accuracy, not stability, is the binding concern, so stay a few x below). Always verify
  with a dt-halving test on <omega^2> (enstrophy) and c_max(t).
- **Resolution requirement:** Delta_x <~ l_a/4. At A = 10000, l_a ~ 0.014, so N = 128
  (Delta_x ~ 0.008) UNDER-resolves (~1.8 pts/l_a) — use N >= 384-512, or lower A to
  ~2000-3000. N = 512 resolves A = 5000-20000 adequately (5-10 pts/l_a).

### Enstrophy diagnostic
`<omega^2>` = mean squared vorticity, omega = -Lap(psi). In Fourier: `omega_tr = k2 .* psi_tr`;
`enstrophy = sum(abs2.(omega_tr))/N^4`. Best convergence diagnostic (dominated by small
scales where errors appear first); pairs with c_max(t) (flow check vs transport check).

---

## 8. The scalar drop initial condition

Current code builds a **tanh flat-top disk** (NOT a Gaussian):
```julia
del = 50/L; r0 = 0.2*L; cx = cy = 0.5*L        # L here = N (grid points)
conc[i,j] = (1 - tanh(del*( sqrt((i-cx)^2+(j-cy)^2) - r0 )))/2
```
In physical units (box length 1): radius r0_phys = 0.2 (diameter = 40% of box),
interface width ~ 2/50 = 0.04. The profile is N-independent (del = 50/L cancels units).

**PROBLEM identified:** r0 = 0.2 is far too big. Conserved mean c_bar ~ pi*r0^2 ~ 0.13,
so c_max can only fall by a factor ~8 (< 1 decade) before the box is uniform — no room
to see any decay law. **FIX: set r0 = 0.05*L** -> c_bar ~ pi*(0.05)^2 ~ 0.008, giving
~2 decades of dynamic range, while keeping r0 ~ 3-5 l_a (multi-vortex regime). Keep the
tanh interface (resolved by ~20 pts at N = 512).

### Always plot c_max - c_bar (mean-subtracted)
In a periodic box c_max decays to the (conserved, analytically known) mean c_bar, not to
zero. **Plot log(c_max - c_bar) vs t.** The late-time slope on a LOG-LINEAR plot gives the
gravest box-mode rate r, from which D_eff = r / (4 pi^2 / L^2) = r L^2 / (4 pi^2).
Without subtraction the floor bends the tail and you measure a spuriously shallow exponent.

### Fitting tip (avoid late-t weighting)
Uniform-in-t sampling over-weights large t on log-log axes. Bin logarithmically in t
before fitting, or weight points by 1/t, or save data at log-spaced times. Always inspect
the *local* slope (d ln C / d ln t) vs t and the compensated plot (c_max - c_bar)*t vs t
(flat where t^-1 holds) rather than trusting a single global fit.

---

## 9. Results obtained so far (A = 5000 assumed for the Diff sweep)

Data files: `C_t__D5.txt`, `C_t__D1.txt`, `C_t__D0_1.txt`, `C_t__D0_01.txt`
(columns: t, c_max - c_bar). Diff = 5, 1, 0.1, 0.01 at fixed A = 5000.

**No clean t^(-1) power-law window exists in these runs.** The local log-log slope sails
through -1 (around t ~ 5-15 for Diff = 5, 1) without forming a plateau and keeps
steepening to -4/-5 — that is the exponential box-mode tail masquerading as an
ever-steepening log-log slope. The box (L = 1, r0 too large) saturates before the
intermediate -1 scaling can develop. **The effective-diffusion PHYSICS is confirmed, but
the box is too small to display the -1 EXPONENT.**

**Effective diffusivity extracted from the late-time exponential (box-mode) rate**
`r_inf`, via `D_eff = r_inf L^2 / (4 pi^2)`:

| Diff | r_inf | D_eff       | D_eff*A | eddy part (D_eff*A - Diff) |
|------|-------|-------------|---------|----------------------------|
| 5    | 0.046 | 1.17e-3     | 5.84    | 0.84 |
| 1    | 0.015 | 3.88e-4     | 1.94    | 0.94 |
| 0.1  | 0.010 | 2.52e-4     | 1.26    | 1.16 |
| 0.01 | (n/a) | (not asymptotic) | -  | -    |

**Predicted vs measured** (model `D_eff = (Diff + const)/A`, single mixing-length
prefactor kappa = 0.45 fixed from the two cleanest tails -> const ~ 0.89):

| Diff | D_eff*A predicted (Diff + 0.89) | D_eff*A measured | ratio |
|------|----------|----------|-------|
| 5    | 5.89     | 5.84     | 0.99  |
| 1    | 1.89     | 1.94     | 1.03  |
| 0.1  | 0.99     | 1.26     | 1.27  |

- Parameter-free molecular check: linear fit D_eff*A = **0.95*Diff + 1.09** vs predicted
  slope 1. The molecular slope (~1) is independently confirmed.
- Velocity ingredient independently verified: predicted u ~ l_a = 0.020 at A = 5000
  matches measured (rescaled) u ~ 0.020.
- Agreement within 1-3% at Diff = 5 and 1 with a single O(1) kappa across a 5x range.

**Honest caveats (do NOT overstate):**
- The eddy part (0.84, 0.94, 1.16) is **NOT flat** — it rises monotonically by ~38% as
  Diff falls. Correct statement: "eddy contribution ~ 0.95, weakly Diff-dependent." The
  upward drift is partly fit bias (low-Diff tails not fully relaxed onto the gravest
  mode -> overestimated rate, worst at small Diff) and partly real (onset of lamellar
  structure as l_B/l_a shrinks; the additive D_mol + D_turb picture is the large-Diff
  limit).
- **Diff = 0.1 oscillates strongly** in its local decay rate (swings 0.01-0.035). This
  is the filamentation crossover: scalar partly stirred into transient filaments
  (sharpen c_max) and partly homogenized (smooth) — neither dominates. Compounded by
  c_max being a single-point (intermittent) statistic and coarser sampling (saved every
  Delta_t = 3). For a clean D_eff from it, fit the smooth lower envelope of the rate.
- **Diff = 0.01 has not reached an asymptotic regime** (c_max still ~0.35 at t = 99) and
  its decay rate is still *increasing* (0.005 -> 0.04) — accelerating decay, the
  signature of approaching the exponential lamellar regime, qualitatively distinct from
  the decelerating approach-to-diffusion of the higher-Diff runs. This is the predicted
  regime change; the run must be extended (t ~ 300-500) to characterize it.

---

## 10. Recommended next steps

1. **Re-run the scalar with r0 = 0.05** (small drop) at N = 512 to open dynamic range,
   and always analyze c_max - c_bar.
2. **A-sweep at fixed Diff = 1** (A = 2500, 5000, 10000, 20000): test that D_turb*A is
   genuinely A-independent (isolates the eddy term; cleaner than the Diff-sweep). This is
   the decisive test of the mixing-length picture and of the "not-flat" eddy trend.
3. **Passive-tracer-cloud measurement of D_eff** via <r^2> = 4 D_eff t (no source, no
   box-mode-fitting ambiguity) as a bias-free cross-check.
4. **Extend Diff = 0.01** (and try Diff = 0.001) to confirm the lamellar regime: decay
   should be exponential in *raw t* at rate ~gamma_dot ~ O(1), NOT following the
   4 pi^2 D_eff / L^2 box-mode law, and curves should collapse in raw t (not t/A).
5. If a clean t^(-1) exponent is required: only A = 20000 with r0 ~ 0.025 opens a ~1-decade
   window (keeps r0 ~ 2.5 l_a). Otherwise lean on the D_eff measurement, which is the
   theory-backed quantity.
6. **Headline figure:** D_eff*A vs Diff — straight line slope 1, intercept ~0.9, with the
   Diff = 0.01 point falling off as the regime changes.
7. Beyond c_max: measure scalar dissipation chi(t) = D<|grad c|^2>, scalar variance,
   the FTLE/stretching-rate distribution, and p(c) (test against Villermaux Gamma /
   convolution predictions) to position active turbulence as a stirring protocol absent
   from the classical mixing catalogue (its stretching statistics are set by defect
   dynamics; control knob gamma_dot ~ S is scale-free; l_a replaces the injection scale).

---

## 11. Quick-reference summary

- 2D active nematic (corotation-only, S = +/-1) + one-way-coupled passive scalar, Stokes flow.
- **A = activity** (large A = more turbulent, finer scales); **R** stabilizes via (1+R/4);
  scalar diffusivity in code is **Diff/A**.
- Active length **l_a = sqrt(2(1+R/4)/A) ~ A^-1/2**; strain rate **gamma_dot ~ S ~ O(1)**.
- **D_eff = (Diff + ~0.9)/A**: molecular + eddy, both ~1/A -> **t/A collapse**. Eddy
  prefactor kappa ~ 0.45. Confirmed to a few % for Diff >= 1; weak residual Diff-dependence.
- Scalar max decays as **t^(-1)** (2D effective diffusion) at high Diff, **exponential**
  (lamellar) at low Diff; box saturation gives exponential tail -> extract D_eff from it.
- Numerics: pseudo-spectral Julia, **N = 512, L = 1**, 2/3 dealias, **IMEX (implicit
  diffusion)**, dt = 0.02, drop r0 = 0.05, plot c_max - c_bar.
- Key reference: **Villermaux, Annu. Rev. Fluid Mech. 2019, 51:245** (lamellar mixing).
