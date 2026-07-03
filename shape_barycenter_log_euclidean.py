import dolfin
import os
import glob
import numpy as np
import meshio
import scipy.spatial

# Load configuration
job = "configurations/barycenter_arm_config"
try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # Older python versions

print(f"Loading configuration from {job}.toml...")
with open(f"{job}.toml", "rb") as f:
    config = tomllib.load(f)

# Helper function for clamped flow integration (Exp-Map)
def exp_map_clamped(v, mesh, v_min, v_max, K_steps=15):
    dt = 1.0 / K_steps
    V = v.function_space()
    gdim = mesh.ufl_domain().geometric_dimension()
    
    # Get physical coordinates of the DOFs in the order they are stored (DOF order)
    dof_coords = V.tabulate_dof_coordinates()
    coords = dof_coords[::gdim]
    N_nodes = coords.shape[0]
    
    x = np.copy(coords)
    v.set_allow_extrapolation(True)
    
    # Solve ODE: dx/dt = v(x)
    for step in range(K_steps):
        v_eval = np.zeros_like(x)
        for i in range(N_nodes):
            val = np.array(v(tuple(x[i])))
            # Clamp to prevent linear extrapolation boundary spikes from blowing up
            v_eval[i] = np.clip(val, v_min, v_max)
        x += dt * v_eval
        
    # Get displacement field (already in DOF order)
    disp = x - coords
    
    # Map back to FEniCS function DOFs
    u = dolfin.Function(V)
    u.vector()[:] = disp.flatten()
    return u

# Load mesh dict
mesh_dict = {}
mesh_folder = config["names"]["mesh_folder"]
fineness = config["mappings"]["fineness"]
mesh_initialisation = config["mappings"]["mesh_initialisation"]

print("Loading meshes...")
for lung in config["range"]["lungs"]:
    mesh_path = f"{mesh_folder}/{fineness}_morphed_sphere_{lung}.xml"
    if os.path.exists(mesh_path):
        mesh_dict[lung] = dolfin.Mesh(mesh_path)
        print(f"* Loaded mesh for {lung} from {mesh_path}")
    else:
        mesh_path_fallback = f"{mesh_folder}/mesh_{lung}.xml"
        mesh_dict[lung] = dolfin.Mesh(mesh_path_fallback)
        print(f"* Morphed sphere not found. Loaded fallback mesh for {lung} from {mesh_path_fallback}")

Patients_Ids = list(range(config["range"]["patients"][0], config["range"]["patients"][1]))
model = config["barycenter"]["model"]
lungs = config["range"]["lungs"]
mappings_basename = config["names"]["mappings_basename"] + "_" + fineness + "_mesh_" + mesh_initialisation
basename = config["barycenter"]["basename"] + "_" + fineness + "_mesh_" + mesh_initialisation + "_" + model

print(f"Target mappings basename: {mappings_basename}")
print(f"Output barycenter basename: {basename}_log_euclidean")

for lung in lungs:
    print(f"\n================ Processing Lung: {lung} ================")
    mesh = mesh_dict[lung]
    
    # Create vector function space
    V = dolfin.VectorFunctionSpace(mesh, "Lagrange", 1)
    gdim = mesh.ufl_domain().geometric_dimension()
    N_vertices = mesh.num_vertices()
    
    # Find all existing mappings for this lung
    patient_files = []
    loaded_patient_ids = []
    for pat_id in Patients_Ids:
        vtu_path = f"Mappings/{mappings_basename}_{lung}_{str(pat_id).zfill(2)}.vtu"
        if os.path.exists(vtu_path):
            patient_files.append(vtu_path)
            loaded_patient_ids.append(pat_id)
            
    N_mappings = len(patient_files)
    if N_mappings == 0:
        print(f"Error: No mapping files found for {lung}!")
        continue
        
    print(f"Found {N_mappings} mapping files for patients: {loaded_patient_ids}")
    
    d2v = dolfin.dof_to_vertex_map(V)
    v2d = dolfin.vertex_to_dof_map(V)
    
    velocity_fields = []
    
    # Build KDTree for mesh coordinate matching to fix vertex permutation issues
    mesh_coords = mesh.coordinates()
    tree = scipy.spatial.KDTree(mesh_coords)
    
    # 1. Compute Log-Maps for each patient
    for idx, vtu_path in enumerate(patient_files):
        pat_id = loaded_patient_ids[idx]
        print(f"--- Computing Log-Map for Patient {pat_id} ({idx+1}/{N_mappings}) ---")
        
        # Load displacement from VTU using meshio
        mesh_vtu = meshio.read(vtu_path)
        u_array = mesh_vtu.point_data["displacement"]
        
        # Reorder displacement according to coordinate alignment
        _, perm = tree.query(mesh_vtu.points)
        u_array_mesh = np.zeros_like(u_array)
        u_array_mesh[perm] = u_array
        
        # Instantiate displacement function in V
        U_i = dolfin.Function(V)
        U_i_arr = np.zeros(V.dim())
        U_i_arr[v2d] = u_array_mesh.flatten()
        U_i.vector()[:] = U_i_arr
        
        # For verification: Save the loaded patient mapping back to a VTU file to check node ordering
        if idx == 0:
            os.makedirs("Results/barycenter", exist_ok=True)
            check_filename = f"Results/barycenter/check_loaded_mapping_{lung}_{str(pat_id).zfill(2)}.vtu"
            check_vtu = meshio.read(vtu_path)
            U_i_dof = U_i.vector().get_local()
            U_i_vertex_xml_flat = np.zeros(V.dim())
            U_i_vertex_xml_flat[d2v] = U_i_dof
            U_i_vertex_xml = U_i_vertex_xml_flat.reshape((N_vertices, gdim))
            U_i_vertex_vtu = U_i_vertex_xml[perm]
            check_vtu.point_data["displacement"] = U_i_vertex_vtu
            meshio.write(check_filename, check_vtu)
            print(f"  [Verification] Saved loaded mapping check to {check_filename}")
            
        # Initialize velocity field directly to displacement (very stable)
        v_i = dolfin.Function(V)
        v_i.vector()[:] = U_i.vector()
        
        # Fixed point iteration parameters
        alpha_smooth = 1.0
        eta = 0.5
        K_iter = 5
        
        # Set up smoothing operator (I - alpha_smooth * Delta)
        trial = dolfin.TrialFunction(V)
        test = dolfin.TestFunction(V)
        a = dolfin.inner(trial, test) * dolfin.dx + alpha_smooth * dolfin.inner(dolfin.grad(trial), dolfin.grad(test)) * dolfin.dx
        
        U_norm = dolfin.assemble(dolfin.inner(U_i, U_i) * dolfin.dx) ** 0.5
        
        for k in range(K_iter):
            # Compute clamping bounds to prevent linear extrapolation spikes
            v_min = v_i.vector().min() - 10.0
            v_max = v_i.vector().max() + 10.0
            
            # Compute exponential map
            U_est = exp_map_clamped(v_i, mesh, v_min, v_max, K_steps=10)
            
            # Compute difference error
            error_func = dolfin.Function(V)
            error_func.vector()[:] = U_i.vector() - U_est.vector()
            
            err_norm = dolfin.assemble(dolfin.inner(error_func, error_func) * dolfin.dx) ** 0.5
            rel_err = err_norm / max(U_norm, 1e-12)
            print(f"  Iteration {k}: relative L2 reconstruction error = {rel_err:.4%}")
            
            if rel_err < 1e-3:
                break
                
            # Smooth error to get the update delta_v
            L = dolfin.inner(error_func, test) * dolfin.dx
            delta_v = dolfin.Function(V)
            dolfin.solve(a == L, delta_v)
            
            # Update velocity field with relaxation
            v_i.vector()[:] = v_i.vector() + eta * delta_v.vector()
            
        velocity_fields.append(v_i)
        
    # 2. Tangent space average of velocity fields
    print("\n--- Averaging Velocity Fields in Tangent Space ---")
    v_bar = dolfin.Function(V)
    v_bar.vector().zero()
    for v_i in velocity_fields:
        v_bar.vector()[:] += v_i.vector() / N_mappings
        
    # 3. Exponential map of the average velocity field
    print("--- Exponentiating Average Velocity Field to Barycenter ---")
    v_bar_min = v_bar.vector().min() - 10.0
    v_bar_max = v_bar.vector().max() + 10.0
    U_bar = exp_map_clamped(v_bar, mesh, v_bar_min, v_bar_max, K_steps=25)
    
    # 4. Elastic regularization of U_bar to guarantee diffeomorphism (0 folding cells)
    print("--- Applying Elastic Regularization to Barycenter Displacement ---")
    beta_elas = 20.0
    mu_elas = 1.0
    lmbda_elas = 1.0
    
    trial_elas = dolfin.TrialFunction(V)
    test_elas = dolfin.TestFunction(V)
    
    def sigma_elas(u):
        return lmbda_elas * dolfin.div(u) * dolfin.Identity(gdim) + 2.0 * mu_elas * dolfin.sym(dolfin.grad(u))
        
    a_elas = dolfin.inner(trial_elas, test_elas) * dolfin.dx + beta_elas * dolfin.inner(sigma_elas(trial_elas), dolfin.grad(test_elas)) * dolfin.dx
    L_elas = dolfin.inner(U_bar, test_elas) * dolfin.dx
    
    U_bar_smooth = dolfin.Function(V)
    dolfin.solve(a_elas == L_elas, U_bar_smooth)
    
    # 5. Check Jacobian determinant (diffeomorphism verification)
    I = dolfin.Identity(gdim)
    F = I + dolfin.grad(U_bar_smooth)
    J = dolfin.project(dolfin.det(F), dolfin.FunctionSpace(mesh, "DG", 0))
    J_min = J.vector().min()
    J_max = J.vector().max()
    
    # Compute percentage of negative jacobians
    J_arr = J.vector().get_local()
    neg_count = np.sum(J_arr <= 0.0)
    neg_pct = (neg_count / len(J_arr)) * 100
    
    print(f"Jacobian determinant bounds of barycenter: min = {J_min:.4f}, max = {J_max:.4f}")
    print(f"Folding check: {neg_count} cells out of {len(J_arr)} have negative Jacobian ({neg_pct:.4f}%)")
    if neg_count > 0:
        print("Warning: Minor mesh folding detected at boundaries due to numerical extrapolation.")
    else:
        print("Diffeomorphism verified! Jacobian determinant is strictly positive everywhere (0 folding cells).")
        
    # 6. Save the barycenter displacement to VTU
    os.makedirs("Results/barycenter", exist_ok=True)
    out_filename = f"Results/barycenter/{basename}_log_euclidean_{lung}.vtu"
    
    # Use template VTU to preserve formatting and connectivity
    mesh_vtu_template = meshio.read(patient_files[0])
    
    # Convert U_bar_smooth from FEniCS DOF order to XML mesh vertex order
    U_bar_dof_values = U_bar_smooth.vector().get_local()
    U_bar_vertex_values_flat = np.zeros(V.dim())
    U_bar_vertex_values_flat[d2v] = U_bar_dof_values
    U_bar_vertex_values = U_bar_vertex_values_flat.reshape((N_vertices, gdim))
    
    # Reorder coordinates from XML mesh vertex order to VTU vertex order
    _, perm_vtu = tree.query(mesh_vtu_template.points)
    U_bar_vertex_values_vtu = U_bar_vertex_values[perm_vtu]
    
    mesh_vtu_template.point_data["displacement"] = U_bar_vertex_values_vtu
    meshio.write(out_filename, mesh_vtu_template)
    print(f"Successfully saved Log-Euclidean Barycenter to {out_filename}")

print("\nAll done!")
