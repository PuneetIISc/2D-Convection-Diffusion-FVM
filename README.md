# 2D Convection–Diffusion Solver (Unstructured, Cell-Centred FVM)

A C++ finite-volume solver for the steady 2D scalar convection–diffusion
equation on **unstructured triangular meshes**, with both **explicit** and
**implicit** time integration. Written as a hands-on study of the finite-volume
method for computational heat- and mass-transfer analysis, and validated
against an exact analytic solution.

<!-- Optional badges — uncomment once you set up the repo / CI
![C++](https://img.shields.io/badge/C%2B%2B-11-blue)
![CMake](https://img.shields.io/badge/build-CMake-informational)
![License: MIT](https://img.shields.io/badge/License-MIT-green)
-->

## Problem

$$
\frac{\partial (u\phi)}{\partial x} + \frac{\partial (v\phi)}{\partial y}
= \alpha\left(\frac{\partial^2\phi}{\partial x^2} + \frac{\partial^2\phi}{\partial y^2}\right)
$$

on the unit square, with velocity $(u,v)=(u_0,0)$ and Péclet number
$Pe = u_0 L/\alpha$.

**Boundary conditions**

| Boundary | Condition |
|---|---|
| Left (inlet) | $\phi = \sin(\pi y)$ |
| Top, bottom (walls) | $\phi = 0$ |
| Right (outlet) | $\partial\phi/\partial x = 0$ |

**Exact steady solution** (used for validation)

$$
\phi^{\text{exact}}(x,y) = \sin(\pi y)
\frac{r_2e^{r_1 x + r_2 L} - r_1e^{r_1 L + r_2 x}}
     {r_2e^{r_2 L} - r_1e^{r_1 L}},
\qquad
r_{1,2} = \frac{u_0}{2\alpha} \pm \sqrt{\frac{u_0^2}{4\alpha^2} + \pi^2}
$$

## Build and run

```bash
cmake -S . -B build
cmake --build build
./build/cd2d userchoice.inp file.input
```

Requirements: a C++11 compiler and CMake ≥ 3.10.

Outputs (written to the working directory): `Solution.dat` (Tecplot),
`Solution.vtk` (ParaView), and `residue.dat` (iteration vs. residue).

## Input files

`file.input` — number of files, then the three mesh files:

```
3
mesh/geometry.out
mesh/nodesupport.out
mesh/cell.out
```

`userchoice.inp` — five lines, read positionally (values first, trailing label ignored):

```
Pe   L   uvel   vvel
order  nodal_value_option  timestep      ! order 1/2 ; nodal 1=area,2=pseudo-Lap ; timestep 1=global,2=local
mode_impli  rkstage  cfl  cfl_option     ! mode 0=explicit RK,1=implicit SGS ; cfl_option 1=const,2=ramp,3=ramp+cap
no_of_sweeps  cflmulti  maxcfl  relax
itermax  maximum_residue
```

## Numerical method

* **Co-volume**: for each face, the diamond area = triangle(n1,n2,incell) + triangle(n1,n2,outcell).
* **Face gradient**: computed by a Green–Gauss contour integral around the co-volume — summing the face-averaged values of $\phi$ times the outward edge normals — and dividing by the co-volume area. This gives $\nabla\phi$ at the face for the diffusion flux.
* **Diffusion**: on each face the diffusive flux is the diffusion coefficient times the gradient projected onto the face's area-normal vector, $\alpha\,\nabla\phi\cdot\mathbf{S}_f$ with $\alpha = 1/Pe$. The gradient is obtained from the Green–Gauss integral over the co-volume, so the scheme stays accurate on skewed cells. At the outlet the normal derivative is set to zero ($\partial\phi/\partial x = 0$), keeping only the tangential variation.
* **Cell gradient**: co-volume-weighted sum of the face diamond gradients.
* **Convection**: first- or second-order upwind reconstruction, upwind flux.
* **Vertex values**: area-weighted or pseudo-Laplacian interpolation from the surrounding cells, with the outlet nodes fixed by a least-squares fit.
* **Time**: explicit multistage Runge–Kutta, or implicit symmetric Gauss–Seidel; local or global time step; $dt = cfl\,\sqrt{\text{volume}}/|\mathbf{v}|$.
* **Residue**: RMS of $|(\phi^{n+1}-\phi^{n})/dt|$, normalised by iteration 1.

## Code layout

```
include/   Vec2, Mesh, Input, Convection, TimeIntegrator, Solver, Writer
src/       Mesh (file reader + co-volumes), Input, Solver,
           TimeIntegrator, Writer, main
mesh/      geometry.out, nodesupport.out, cell.out
```

The mesh is stored **face-by-face** (each `Face` knows its two nodes, its
in/out cells, and its area-normal), so the solver is a single loop over faces —
the same organisation used by OpenFOAM and SU2. Convection schemes and time
integrators are small class hierarchies, so adding one is just a new derived
class.

## Validation

All runs use $u_0 = 1$, $L = 1$ (so $\alpha = 1/Pe$) and are compared against the
exact solution through the area-weighted global $L_2$ error.

### 1. Solution field (Pe = 100)

The computed field reproduces the exact solution: a sinusoidal profile carried
across the domain and thinned toward the outlet. On the coarse mesh the contours
are faceted; on the fine mesh the numerical and exact fields are visually
identical.

<!-- Add images to a docs/ folder and they will render here:
![Solution, fine mesh](docs/fine_sol.png)
-->

### 2. Grid convergence (order of accuracy)

The global $L_2$ error was measured on four successively refined meshes
(Pe = 100, second-order upwind, pseudo-Laplacian nodes):

| Points (np) | Cells (nc) | $L_2$ error |
|---:|---:|---:|
| 68   | 106  | 4.902 × 10⁻³ |
| 241  | 424  | 1.668 × 10⁻³ |
| 905  | 1696 | 6.049 × 10⁻⁴ |
| 3505 | 6784 | 2.026 × 10⁻⁴ |

A least-squares fit of $\log(L_2)$ vs. $\log(h)$, with $h = 1/\sqrt{n_c}$, gives an
**observed order of ≈ 1.53**. This lies between first and second order, which is
the documented behaviour of this viscous discretisation on distorted triangles
(locally inconsistent but globally convergent, see Ref. [1]); against the
reference measure $dx = 1/n_p$ the same data reproduce the fall rate reported
there.

<!-- ![Grid convergence](docs/order_of_convergence.png) -->

### 3. Effect of the Péclet number

Run on the finest mesh at Pe = 1, 10 and 100. As Pe increases the contours
change from rounded and symmetric (diffusion-dominated) to stretched in the
flow direction with a thin outlet layer (convection-dominated) — captured
cleanly, with no spurious oscillations.

<!-- ![Peclet study](docs/isoplot.png) -->

### 4. Explicit vs. implicit

Same spatial scheme, so both reach the **same steady state**; only the
convergence speed differs:

| Scheme | Iterations to `1e-10` |
|---|---:|
| Explicit (3-stage Runge–Kutta, CFL 0.2) | 190 |
| Implicit (symmetric Gauss–Seidel, ramped CFL) | 23 |

The residual falls monotonically in both cases, and a line plot of $\phi$ through
the domain shows the explicit and implicit curves lying exactly on top of one
another — an ≈ 8× reduction in iterations at no cost in accuracy.

<!-- ![Residual convergence](docs/residual_convergence.png) -->

> A full write-up of the method, cases and results is in
> [`docs/report/CD2D_Technical_Report.pdf`](docs/report/CD2D_Technical_Report.pdf).

## References

1. N. Munikrishna, *On Viscous Flux Discretization Procedures for Finite Volume
   and Meshless Solvers*, Ph.D. Thesis, Department of Aerospace Engineering,
   Indian Institute of Science, Bengaluru.
2. F. Moukalled, L. Mangani and M. Darwish, *The Finite Volume Method in
   Computational Fluid Dynamics: An Advanced Introduction with OpenFOAM and
   Matlab*, Springer, 2016.

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
