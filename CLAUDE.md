# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

This is a Julia-based active-nematic fluid dynamics solver for a TFG (undergraduate thesis). The solver simulates **active turbulence** in 2D using a pseudo-spectral (Fourier-space) method with explicit Euler time integration.

The primary file is `STUDENTS_defect_free_turbu.ipynb` — a Jupyter notebook running Julia via IJulia.

## Running the notebook

Open and run cells in Jupyter with a Julia kernel (IJulia). There is no build step or test suite. The notebook is the entire codebase.

To restart a simulation from a saved state, uncomment and adjust the loading block:
```julia
theta = reshape( readdlm("results/julia/L_512_A_1000000/last_theta_101.txt"), (L,L) )
```

## Physics and numerical architecture

The solver evolves the **nematic director angle** `theta(x,y,t)` on a 2D periodic domain of size `L×L`.

**Two-function structure:**
- `diff_step(conc_tr)` — simple diffusion test (Fourier-space: multiply by `-k²`)
- `turbu_step(theta, A, R, S)` — full active-nematic step; returns `(theta_out, theta_tr, psi_tr)`

**`turbu_step` pipeline:**
1. Dealias and FFT `theta`
2. Compute `sin(2θ)`, `cos(2θ)` (nonlinear terms, done in real space)
3. Build the Ericksen stress `eric` (cross-product of gradient and biharmonic gradient)
4. Solve for stream function `psi` in Fourier space by inverting `∇⁴` (regularized as `k4_im` to avoid division by zero at k=0)
5. Compute covariant derivative of `theta` via `psi`
6. Euler step on `theta_tr`

**Key parameters:**
- `A` — activity (large A → turbulent regime; tested at `A=100000`)
- `R` — viscosity ratio
- `S` — `+1` contractile, `-1` extensile

**Fourier operators** (all global, set once before the time loop):
- `kx`, `ky` — wavenumber arrays (shaped for broadcasting: `(L,1,1)` and `(1,L,1)`)
- `k2`, `k4`, `k4_im` — squared and quartic wavenumber magnitudes
- `dealias` — 2/3 dealiasing mask (zeros out modes with `|k| > k_max`)

**Output layout:** results saved to `L_<L>_A_<A>/` with subdirectories `theta/` and `U/`, plus `params.txt`. The last theta frame is also saved as `last_theta_<n>.txt` for restarts.

## Context management

`MAINTAIN_prompt.md` contains a recurring prompt to keep a `CONTEXT.md` file in sync with the notebook at the end of each working session. If a `CONTEXT.md` exists in the repo, treat it as the authoritative living document for current solver state, known issues, and next steps — it takes precedence over notebook comments for high-level context.
