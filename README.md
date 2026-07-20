# Physics Simulations Portfolio

Six C++/Python projects applying numerical methods to physics and quantitative-finance problems:
root-finding, explicit and symplectic time integration, finite-difference PDE solvers, and rigid-body
mechanics.

Each project is self-contained: a short, focused C++ program (or a few), a Jupyter notebook that reads
its CSV output and produces the plots below, and a `README.md` explaining the problem and the method.

## Projects

| Project | What it does |
|---|---|
| [`finite-difference-pde/`](finite-difference-pde/) | Explicit finite-difference solver for the Black-Scholes PDE on a non-uniform grid (call + put, with a full stability analysis), plus a 1D heat-equation solver |
| [`n-body-solar-system/`](n-body-solar-system/) | Velocity-Verlet N-body gravity simulation, from two bodies to a nine-body solar system, verified energy-conserving to ~1 part in 10¹⁰ |
| [`irregular-body-collision-2d/`](irregular-body-collision-2d/) | A small 2D rigid-body physics engine: an arbitrary shape rasterized from a grid, falling and bouncing with a closed-form elastic/inelastic collision response that exchanges translational and rotational energy |
| [`sphere-collisions/`](sphere-collisions/) | Bisection root-finding applied to elastic 2D collisions, from head-on equal masses up to an oblique-angle pool-table scenario |
| [`bar-linkage-kinematics/`](bar-linkage-kinematics/) | Closed-form position analysis and finite-difference velocity/acceleration of a 5-bar planar linkage |
| [`projectile-motion/`](projectile-motion/) | Numerical-vs-analytic convergence analysis, linear-drag trajectories, and a two-projectile interception search |

## Why these projects

Every project here is really the same underlying skill applied to a different equation: turn a
continuous physical model into a discrete numerical scheme, get the scheme's correctness conditions
right (stability limits, valid root-finding brackets, symplectic integration), and validate the result
against a conservation law or a closed-form limit rather than just "it looks right." That's the same
mindset that error analysis, PDE-based pricing, and time-series simulation call for in quantitative
finance — `finite-difference-pde/black-scholes` is the most direct example, but the validation habit
(checking energy conservation, checking convergence order, deriving valid parameter ranges instead of
guessing them) shows up in every project.

## Notes

- Each project's `README.md` has the exact build command; all C++ code targets `g++ -std=c++17`, no
  external dependencies beyond the C++ standard library. Notebooks use pandas/matplotlib/numpy.
- Licensed under MIT (see `LICENSE`).
