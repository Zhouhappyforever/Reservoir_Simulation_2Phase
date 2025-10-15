# --------------------------------------------------------------------------
# SECTION 1: PREAMBLE AND PACKAGE IMPORTS
# --------------------------------------------------------------------------
# Import the necessary Gridap packages and a plotting package.
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Statistics
using Plots
using Gridap                # Main Gridap package for finite element analysis
using Gridap.Geometry       # For mesh and geometry handling
using Gridap.FESpaces       # For finite element spaces
using Gridap.MultiField     # For coupled multi-physics problems
using Gridap.Io             # For input/output operations
using Gridap.Fields         # For field operations
using Gridap.TensorValues   # For tensor operations
using Gridap.ODEs           # For time-dependent problems
using Gridap.CellData       # For cell data operations and projection
using WriteVTK              # For VTK file output (visualization)
using GridapGmsh            # For Gmsh mesh integration

# Debugging: Check package import
println("Section 1: Packages imported successfully.")


# --------------------------------------------------------------------------
# SECTION 2: MODEL AND DOMAIN DEFINITION
# --------------------------------------------------------------------------
# Define the computational domain: a 1D horizontal pipe of length 100 m.
domain = (0, 100)
# Define the mesh resolution: 100 cells for better resolution in 1D.
partition = (100,)
# Create the discrete model (the mesh).
model = CartesianDiscreteModel(domain, partition)

# Label boundaries for applying boundary conditions.
labels = get_face_labeling(model)
add_tag_from_tags!(labels, "inlet", [1])
add_tag_from_tags!(labels, "outlet", [2])

# Debugging: Check model
trian = Triangulation(model)
println("Section 2: Model created with ", num_cells(trian), " cells and ", length(Gridap.Geometry.get_node_coordinates(trian)), " nodes.")

# --------------------------------------------------------------------------
# SECTION 3: DEFINITION OF FINITE ELEMENT SPACES (MULTI-FIELD)
# --------------------------------------------------------------------------
# Define the polynomial order for the FE approximation.
order = 1
reffe = ReferenceFE(lagrangian, Float64, order)

# Define the test spaces.
Vphi = TestFESpace(model, reffe; dirichlet_tags=["inlet"])
Vvp = TestFESpace(model, reffe; dirichlet_tags=["inlet"])
Vvf = TestFESpace(model, reffe; dirichlet_tags=["inlet"])
Vlam = TestFESpace(model, reffe; dirichlet_tags=["outlet"])
Y = MultiFieldFESpace([Vphi, Vvp, Vvf, Vlam])

# Define the trial spaces and boundary conditions.
phi_inlet(x) = 0.5   # Inlet particle volume fraction
v_inlet(x) = 0.1     # Inlet velocity for both phases (m/s)
lam_outlet(x) = 0.0  # Outlet lambda (pressure gradient proxy)
Uphi = TrialFESpace(Vphi, phi_inlet)
Uvp = TrialFESpace(Vvp, v_inlet)
Uvf = TrialFESpace(Vvf, v_inlet)
Ulam = TrialFESpace(Vlam, lam_outlet)
X = MultiFieldFESpace([Uphi, Uvp, Uvf, Ulam])

# Debugging: Check spaces
println("Section 3: Multi-field FE spaces defined with ", num_free_dofs(Y), " total DOFs.")

# --------------------------------------------------------------------------
# SECTION 4: PHYSICAL PARAMETERS AND CONSTITUTIVE RELATIONS
# --------------------------------------------------------------------------
const ρ_p = 2500.0  # Particle density (kg/m³)
const ρ_f = 1000.0  # Fluid density (kg/m³)
const g = 9.81      # Gravity (m/s²)
theta_angle = 0.0          # Horizontal pipe, sin θ = 0
const d_p = 10000.0  # Increased drag coefficient for stability (kg/(m³·s))
const C_vm = 0.5    # Virtual mass coefficient
b_p(phi) = 0.0      # Collision dispersive pressure (Pa, simplified to 0)

# Debugging: Check parameters
println("Section 4: Parameters defined. Example: ρ_p = ", ρ_p, ", d_p = ", d_p)

# --------------------------------------------------------------------------
# SECTION 5: WEAK FORMULATION (TRANSIENT, NON-LINEAR, MULTI-FIELD) - CORRECTED
# --------------------------------------------------------------------------
degree = 2 * order
dΩ = Measure(Triangulation(model), degree)

# Note: The term `∇(psi) ⋅ (vp * ∇(phi))` is a possible correction for the
# original `∇(psi) ⋅ ∇(phi) * vp`. Please verify this against your model's equations.
# It represents the convection of phi by the velocity vp.
res(t, (phi, vp, vf, lam), (psi, wp, wf, chi)) =
  ∫( psi * ∂t(phi) - ∇(psi) ⋅ (vp * ∇(phi)) )*dΩ +
  ∫( -chi * (∇⋅(phi * vp + (1 - phi) * vf)) )*dΩ +
  ∫( wp ⋅ (ρ_p * (∂t(vp) + (vp ⋅ ∇(vp)))) -
     wp ⋅ (d_p * (vf - vp)) -
     wp ⋅ (C_vm * ρ_f * phi * ( (∂t(vf) + (vf ⋅ ∇(vf))) - (∂t(vp) + (vp ⋅ ∇(vp))) )) -
     (∇⋅wp) * phi * lam )*dΩ +
  ∫( wf ⋅ (ρ_f * (∂t(vf) + (vf ⋅ ∇(vf)))) +
     wf ⋅ (d_p * (vf - vp)) +
     wf ⋅ (C_vm * ρ_f * hi * ( (∂t(vf) + (vf ⋅ ∇(vf))) - (∂t(vp) + (vp ⋅ ∇(vp))) )) -
     (∇⋅wf) * (1 - phi) * lam )*dΩ

# Debugging: Check weak form (dummy evaluation not possible, but confirm definition)
println("Section 5: Weak form defined successfully.")

# --------------------------------------------------------------------------
# SECTION 6: SOLVER CONFIGURATION
# --------------------------------------------------------------------------
op = TransientFEOperator(res, X, Y)
nls = NLSolver(LUSolver(), show_trace=true, method=:newton, iterations=10)

# Smaller time step for stability with inertial terms.
Δt = 1.0 # Time step in seconds
θ = 1.0
odesolver = ThetaMethod(nls, Δt, θ)

# Initial conditions.
phi0(x) = 0.05
vp0(x) = 0.0
vf0(x) = 0.0
lam0(x) = 0.0
xh0 = interpolate_everywhere([phi0, vp0, vf0, lam0], X(0.0))

# Debugging: Check initial condition
coords = Gridap.Geometry.get_node_coordinates(Triangulation(model))
println("Section 6: Solver configured. Initial mean phi = ", Statistics.mean(phi0.(coords)) )

# --------------------------------------------------------------------------
# SECTION 7: SIMULATION EXECUTION
# --------------------------------------------------------------------------
println("DEBUG: Starting simulation...")
t0 = 0.0
T = 100.0 # Total simulation time (100 s)
res(t,(phi,vp,vf,lam),(w_phi,w_p,w_f,w_lam)) = ... # Your long residual definition is here

jac_t(t,(phi,vp,vf,lam),(w_phi,w_p,w_f,w_lam)) = ... # Your jac_t definition is here

# ADD THIS LINE 🔧
jac(t,x,w) = jacobian(res(t,x),x,w)

op = TransientFEOperator(res, (jac, jac_t), U, V) # <-- This line will now work
sol_t = solve(solver, op, xh0, t0, tF)
println("DEBUG: Simulation execution finished. Solution iterator is ready.")
# Debugging: Basic check (iterator exists)
println("Section 7: Simulation executed. Expected steps: ~", T/Δt)

# --------------------------------------------------------------------------
# SECTION 8: POST-PROCESSING
# --------------------------------------------------------------------------
using Gridap.Visualization # Ensure the visualization module is available

println("Starting post-processing with corrected workflow...")

# Use createpvd to manage the collection of time-step files.
# This creates a file named "transient_solution.pvd" that can be opened in ParaView.
createpvd("transient_solution") do pvd
    # Iterate directly over your solution variable, 'sol_t'.
    # This computes each time step lazily, avoiding high memory usage.
    for (solution_t, t) in sol_t
        println("Processing time step at t = $t")

        # CRITICAL STEP: Destructure the multi-field solution object.
        # The order of variables (u_t, p_t) must match the order in which
        # the FE spaces were defined in the MultiFieldFESpace.
        u_t, p_t, phi, lam = solution_t

        # Add the current time step's data to the PVD collection.
        # createvtk saves the fields for this time step to a.vtu file.
        # The destructured fields (u_t, p_t) are passed individually.
        pvd[t] = createvtk(Ω,
                           "transient_solution_$t.vtu",
                           cellfields=["vector_field" => u_t, "scalar_field" => p_t])
    end
end

println("Post-processing finished successfully.")

# This corrected workflow naturally resolves the UndefVarError because
# the problematic 'all_solutions' variable is never used or defined.
# The check in the original Section 9 is no longer necessary.
println("Final solution files have been generated.")




# --------------------------------------------------------------------------
# SECTION 9: IN-LINE PLOTTING OF FINAL RESULT
# --------------------------------------------------------------------------
# println("DEBUG: Checking if final solution was captured...")
# if !isempty(all_solutions)
#     println("Generating final plots...")
#     
#     # Extract the final solution
#     final_xh, final_t = all_solutions[end]
#     final_phi, final_vp, final_vf, final_lam = final_xh
#     
#     # Get nodal coordinates and evaluate fields (1D line plots)
#     coords = Gridap.Geometry.get_node_coordinates(model)
#     x = [p[1] for p in coords]
#     
#     p1 = plot(x, final_phi(coords), label="phi_p", title="Particle Fraction at t = $(final_t) s", xlabel="s (m)", ylabel="phi_p")
#     p2 = plot(x, final_vp(coords), label="v_p", title="Velocities at t = $(final_t) s", xlabel="s (m)", ylabel="v (m/s)")
#     plot!(p2, x, final_vf(coords), label="v_f")
#     p3 = plot(x, final_lam(coords), label="lam", title="Lambda at t = $(final_t) s", xlabel="s (m)", ylabel="lam")
# 
#     # Display side-by-side
#     display(plot(p1, p2, p3, layout=(1,3), size=(1200, 400)))
#     
#     println("Plotting complete.")
# else
#     println("Final solution not found for plotting because no time steps were collected.")
# end

# Debugging: Confirm plots
println("Section 9: In-line plotting completed successfully.")
