# 2D Channel Flow (Plane Poiseuille) – Finite Volume Solver in MATLAB


**Author:** Rochitha A T (AM24S012)

A MATLAB implementation of a 2D, steady, incompressible, laminar channel flow (Re = 100) solved with a Finite Volume Method on a structured staggered grid, with PISO pressure–velocity coupling. The results are compared against the analytical Plane Poiseuille solution.

---

## Problem Definition

A Newtonian fluid enters a 2D channel of height `H = 1` and length `L = 10` with a uniform velocity `u = 1` (Re = 100). The task is to solve the 2D incompressible Navier–Stokes equations and:

1. Plot the streamlines
2. Plot the velocity magnitude
3. Plot the pressure distribution
4. Plot the u-velocity profiles at x = 1, 2, …, 10
5. Find the entry length `L_E` where the flow becomes fully developed (du/dx = 0)
6. Compare with the analytical Plane Poiseuille solution

### Governing equations

Continuity:

```
∂u/∂x + ∂v/∂y = 0
```

Momentum (steady, laminar, constant viscosity):

```
ρ (u ∂u/∂x + v ∂u/∂y) = -∂p/∂x + μ (∂²u/∂x² + ∂²u/∂y²)
ρ (u ∂v/∂x + v ∂v/∂y) = -∂p/∂y + μ (∂²v/∂x² + ∂²v/∂y²)
```

For fully developed flow this reduces to `μ d²u/dy² = dp/dx`.

### Boundary conditions

| Boundary | Condition |
|---|---|
| Inlet (x = 0) | u = 1, v = 0, ∂p/∂x = 0 |
| Outlet (x = 10) | ∂u/∂x = 0, ∂v/∂x = 0, p = p_atm |
| Top wall (y = 1) | u = 0, v = 0, ∂p/∂x = 0 |
| Bottom wall (y = 0) | u = 0, v = 0, ∂p/∂x = 0 |

---

## Numerical Method

- **Discretisation:** Finite Volume Method on a structured grid (cell-centred pressure and velocity, face velocities stored separately)
- **Diffusion terms:** second-order central differencing
- **Convection terms:** central differencing in the current implementation
- **Time marching:** explicit, with time step `dt`
- **Pressure–velocity coupling:** PISO
  - Momentum predictor (H-operator, pressure gradient excluded)
  - Two corrector steps (`n_correctors = 2`)
  - Pressure-correction equation solved with SOR (`omega = 1.5`, up to 50 iterations)
- **No under-relaxation**
- **Convergence:** max change of u, v, p between iterations below `tol`

Default numerical parameters: `Nx = 100`, `Ny = 20`, `dt = 0.01`, `max_iter = 10000`, `tol = 1e-6`.

---



---

## Input Files

### `AM24S012_inputfile.xlsx` (solver)

Sheet with columns `Problem Parameters` and `Value` (range `A1:C10`). Required rows:

| Problem Parameters | Example Value |
|---|---|
| Channel length | 10 |
| Channel height | 1 |
| Inlet velocity | 1 |
| Reynolds number | 100 |
| Density | 1 |
| Number of cells in x-direction | 100 |
| Number of cells in y-direction | 20 |
| Time step | 0.01 |
| Maximum number of iterations | 10000 |
| Convergence tolerance | 1e-6 |

### `AM24S012_streamline_input.xlsx` (streamlines)

Columns `parameter` and `value` with rows: `Channel height`, `Channel length`, `Reynolds number`, `Number of points in x-direction`, `Number of points in y-direction`.

---

## How to Run

**Requirements:** MATLAB (R2019b or later recommended). No additional toolboxes are needed beyond base MATLAB.

1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/<your-repo>.git
   cd <your-repo>
   ```
2. Open MATLAB in the repository folder and make sure the `.xlsx` input files are in the same directory as the scripts.
3. Edit the input files if you want different parameters.
4. Run:
   ```matlab
   channel_flow_solver      % flow solution, profiles, entry length, analytical comparison
   streamline_plot          % streamlines
   mesh_plot                % structured mesh
   ```

---

## Post-processing Outputs

- Velocity magnitude contour (numerical vs. analytical)
- Pressure distribution contour (numerical vs. analytical)
- u-velocity profiles at x = 1…10, compared with the analytical parabolic profile `u = 6 u_in y (H − y) / H²`
- Average |du/dx| along the channel and the detected entry length
- Streamline plot with inset showing inlet, developing and fully developed profiles
- Structured mesh plot

---

## Results

| Quantity | Value |
|---|---|
| Entry length `L_E` (numerical, du/dx criterion) | 9.85 units (`L_E/H` = 9.85) |
| Theoretical entry length (0.06 · Re · H) | 6.00 units |

Add your figures to `results/` and reference them here, for example:

```markdown
![Velocity magnitude](results/velocity_magnitude.png)
![Pressure distribution](results/pressure.png)
![U-velocity profiles](results/u_profiles.png)
![Streamlines](results/streamlines.png)
![Mesh](results/mesh.png)
```

---

## Known Limitations

- The streamline script builds a **prescribed** velocity field (linear blend from the uniform inlet profile to a parabola over `0.05·Re·H`) and integrates it into a stream function. It does not use the solved Navier–Stokes field.
- The numerical solution does not yet reach the analytical fully developed profile (peak velocity of 1.5 · u_in) within the domain, and the numerical pressure magnitude differs from the analytical linear pressure drop. The numerical entry length (9.85) is also larger than the theoretical estimate (6.0). Possible causes include the explicit time-stepping, the pressure boundary treatment and the cell-centre velocity update; these are areas for improvement.
- The "analytical" velocity/pressure fields in the entry region use a simple linear blend, not an exact developing-flow solution.

---

## Possible Improvements

- First-order upwind or higher-order (QUICK/TVD) convection schemes
- Proper Rhie–Chow interpolation and consistent pressure boundary conditions
- Mass-conserving outlet treatment
- Mesh refinement and grid-independence study
- Use the solved velocity field to compute the true streamfunction

---

## References

- White, F. M., *Fluid Mechanics* – Plane Poiseuille flow
- Versteeg, H. K. and Malalasekera, W., *An Introduction to Computational Fluid Dynamics: The Finite Volume Method*
- Issa, R. I., "Solution of the implicitly discretised fluid flow equations by operator-splitting", *J. Comput. Phys.*, 1986 (PISO)

---

