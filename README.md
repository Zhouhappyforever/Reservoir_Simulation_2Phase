<<<<<<< HEAD
# PGE 383 Advanced Geomechanics Final Project 

## Wellbore Poroelasticity Simulation

### Overview
In this final project, you will implement a poroelastic simulation of a wellbore in plane strain conditions. You will use the [Gridap.jl](https://github.com/gridap/Gridap.jl) finite element framework to solve the coupled equations of poroelasticity around a wellbore geometry.

You will simulate the poroelastic response of rock around a wellbore. The wellbore is subjected to a pressure differece from the initial pore pressure in the formation. This creates a pressure gradient that drives fluid flow and causes deformation of the rock matrix.

### Assignment Tasks

1. Complete the implementation of the `wellbore.jl` file using `mandel.jl` as a reference.
    - You can run the `mandel.jl` file with the terminal command `juila --project mandel.jl`
2. Use the provided `wellbore.geo` file as the input geometry.
3. Apply the following boundary conditions:
   - At the `top_bottom` boundaries: Fix both $x$ and $y$ displacements ($u_x = 0, u_y = 0$)
   - At the `wellbore` boundary: Apply a constant pressure of $p_b = 30.5$ MPa
   - Initial pressure throughout the domain: $p_0 = 20.0$ MPa
4. Use the material properties and simulation parameters already defined in `wellbore.jl`.
5. Implement the weak form of the poroelastic equations.
6. Set up the solver and time-stepping scheme.
7. Generate output files for visualization.

### Submission Requirements
- Submit your completed `wellbore.jl` file via Canvas.
- Your code should run without errors when executed with the provided `wellbore.geo` file.

### Evaluation Criteria
- Correctness of the implementation
- Proper application of boundary conditions
- Code organization and clarity
- Proper implementation of the weak form
- Successful execution of the simulation

### Note
All submitted code will be executed to verify functionality. Make sure your code runs correctly with the provided geometry file and parameters.

### Visualization in ParaView

You can use the open source tool [Paraview](https://www.paraview.org/) for visualization.
Unfortunately, you cannot use this from our web-based development environment so
you must download the files to you local computer.  First, on the Terminal
command line in the repository directory, run

```bash
zip -r results* ./results/*
```

Once you do that, in the file browser, right-click `results.zip` and select
"Download".  This will download the file your local computer.  Once you have it
downloaded, you can unzip it locally, either by running

```bash
unzip results.zip
```

in a Terminal session or by simply double-clicking the file on most operating
systems.  Now you're ready to open them in Paraview.

If you don't have Paraview installed, it can be downloaded [here](https://www.paraview.org/download/) for all major operating systems.

Open the `results.pvtu` file from within Paraview and click "Apply" in the Properties tab on the left-hand side of the screen.  Then you should be able to change the selection from "Solid Color" to "pressure" or "displacement".
=======
# Reservoir_Simulation_2Phase
The weak form of the governing equations is derived for the 1D wellbore flow model in Frenet coordinates, assuming a curved wellbore parametrized by arc length \(s \in [0, L]\). The system consists of mass balances for the particle and fluid phases (with volume fraction constraint \(\phi_f = 1 - \phi_p\)) and momentum balances for each phase. The weak form is obtained by multiplying each strong equation by a test function, integrating over the domain, and applying integration by parts where appropriate to reduce the order of derivatives (suitable for finite element discretization). Boundary terms are assumed to vanish for simplicity (e.g., due to no-flux or periodic boundary conditions; in practice, these would be incorporated based on specific wellbore inlet/outlet conditions).

We seek solutions \(\phi_p(s, t)\), \(v_p(s, t)\), \(v_f(s, t)\), and \(\lambda(s, t)\) in appropriate function spaces (e.g., Sobolev spaces \(H^1([0, L])\) for variables with gradients), such that for all test functions \(\psi(s)\), \(\chi(s)\), \(w_p(s)\), and \(w_f(s)\) in the corresponding test spaces (e.g., \(H^1_0([0, L])\) to handle vanishing boundaries), the following variational equations hold:

### Weak Form of Mass Balance for Particle Phase
\[
\int_0^L \psi \frac{\partial \phi_p}{\partial t} \, ds - \int_0^L \frac{\partial \psi}{\partial s} (\phi_p v_p) \, ds = 0.
\]

### Weak Form of Mass Balance for Fluid Phase
\[
\int_0^L \chi \frac{\partial \phi_f}{\partial t} \, ds - \int_0^L \frac{\partial \chi}{\partial s} (\phi_f v_f) \, ds = 0,
\]
where \(\phi_f = 1 - \phi_p\). (Note: This equation is dependent on the particle mass balance due to the saturation constraint, but both are included for completeness. In numerical implementations, one may be eliminated, and the constraint enforced pointwise.)

### Weak Form of Momentum Balance for Particle Phase
\[
\int_0^L w_p \rho_p \left( \frac{\partial v_p}{\partial t} + v_p \frac{\partial v_p}{\partial s} \right) \, ds - \int_0^L w_p \rho_p g \sin \theta(s) \, ds - \int_0^L w_p d_p (v_f - v_p) \, ds \\
- \int_0^L w_p C_{vm} \rho_f \phi_p \left( \frac{\partial v_f}{\partial t} + v_f \frac{\partial v_f}{\partial s} - \frac{\partial v_p}{\partial t} - v_p \frac{\partial v_p}{\partial s} \right) \, ds \\
+ \int_0^L \frac{\partial w_p}{\partial s} \phi_p \lambda \, ds + \int_0^L \frac{\partial w_p}{\partial s} b_p \, ds = 0,
\]
where \(b_p = b_p(\phi_p)\) is the collision dispersive pressure (e.g., \(b_p \propto \phi_p^2 a^2 (\dot{\gamma})^2\), with \(\dot{\gamma}\) the local shear rate).

### Weak Form of Momentum Balance for Fluid Phase
\[
\int_0^L w_f \rho_f \left( \frac{\partial v_f}{\partial t} + v_f \frac{\partial v_f}{\partial s} \right) \, ds - \int_0^L w_f \rho_f g \sin \theta(s) \, ds + \int_0^L w_f d_p (v_f - v_p) \, ds \\
+ \int_0^L w_f C_{vm} \rho_f \phi_p \left( \frac{\partial v_f}{\partial t} + v_f \frac{\partial v_f}{\partial s} - \frac{\partial v_p}{\partial t} - v_p \frac{\partial v_p}{\partial s} \right) \, ds \\
+ \int_0^L \frac{\partial w_f}{\partial s} \phi_f \lambda \, ds = 0.
\]

### Notes on Derivation and Implementation
- **Integration by Parts**: Applied to the advective terms in the mass balances (to shift derivatives from the flux to the test function) and to the pressure (\(\partial \lambda / \partial s\)) and dispersion (\(\partial b_p / \partial s\)) terms in the momentum balances (to shift derivatives to the test functions, facilitating \(C^0\) finite element approximations).
- **Physical Interpretation**: The weak forms preserve the conservation properties of the strong equations while allowing for discontinuous solutions (e.g., shocks in hyperbolic regimes) and numerical stabilization (e.g., SUPG for advection-dominated terms).
- **Coupling and Constraints**: The volume fraction constraint \(\phi_p + \phi_f = 1\) is enforced pointwise. The Lagrange multiplier \(\lambda\) (mixture pressure) couples the momentum equations; the combined pressure contributions across phases yield \(\int_0^L \lambda \left( \phi_p \frac{\partial w_p}{\partial s} + \phi_f \frac{\partial w_f}{\partial s} \right) ds\), reflecting the incompressibility condition \(\phi_p \frac{\partial v_p}{\partial s} + \phi_f \frac{\partial v_f}{\partial s} = 0\) in the steady limit.
- **Limitations**: The forms assume smooth \(\theta(s)\) (local inclination); curvature \(\kappa(s)\) and torsion \(\tau(s)\) do not appear explicitly but may require additional terms for secondary flows in highly curved sections. Time discretization (e.g., implicit Euler) would treat the \(\partial / \partial t\) terms as mass matrix contributions in finite elements.
- **Numerical Considerations**: Suitable for Galerkin finite elements with linear or quadratic basis functions. Stabilization may be needed for high Peclet numbers (advection vs. dispersion) or Reynolds numbers (inertia vs. drag). Boundary conditions (e.g., inlet velocity, outlet pressure) would modify the ignored boundary terms.
>>>>>>> d72ada8a07f6064467d25dfbad8f2a81e1f5d5bb
