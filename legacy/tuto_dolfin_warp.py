import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import dolfin
import dolfin_warp as dwarp
import mshr
import dolfin_mech as dmech
# import create_data
import glob


#%% Domain omega generation

sphere_center   = (120, 120, 120)                                                           # Center of the sphere
sphere_radius   = 75                                                                        # Radius of the sphere
resolution      = 30                                                                        # Resolution of the mesh

# Create a 3D spherical domain
center          = dolfin.Point(sphere_center[0], sphere_center[1], sphere_center[2])        # Center of the disc
radius          = 1.4*sphere_radius                                                         # Radius of the disc
domain          = mshr.Sphere(center, radius)
mesh_omega      = mshr.generate_mesh(domain, resolution)
mesh_omega.num_vertices()

#%% load previous lung mesh


mesh_lung=dolfin.Mesh("Meshes/mesh_RL.xml")


#%% Load image

mesh_name                   = "3D_lung_PA5"
image_basename              = "PA5_Binary"
image_suffix                = "signed_int"
result_folder               = "Results" 
filebasename                = "Lagrange_fine_Grad_Backtracking"
image_name                  = image_basename+"_"+image_suffix
image_folder                = "Images"
reduced_kinematics_model    = "translation+scaling"
image_path                  = image_folder+"/"+image_name


# Checks folders existence
import os 
if not os.path.isdir(image_folder):
    os.system("mkdir "+image_folder)
    print("folder "+image_folder+" created")

if not os.path.isdir(result_folder):
    os.system("mkdir "+result_folder)
    print("folder "+result_folder+" created")

image_file = glob.glob(image_path)


#%% shape registration via dwarp
import time
t_start = time.time()
dwarp.warp(
        warping_type                                = "registration"                                            , # registration
        working_folder                              = result_folder,
        images_char_func                            = False,
        working_basename                            = filebasename,
        images_folder                               = image_folder,
        images_basename                             = image_name,
        images_ext                                  = "vti",
        mesh                                        = mesh_omega,
        kinematics_type                             = "full",
        reduced_kinematics_model                    = reduced_kinematics_model,
        images_quadrature                           = 1,
        n_iter_max                                  = 15,
        nonlinearsolver                             = "gradient_descent",
        min_gradient_step                           = 1e-6,
        gradient_step                               = 1,
        continue_after_fail                         = 1,
        write_VTU_files                             = True,
        write_VTU_files_with_preserved_connectivity = True,
        initialize_reduced_U_from_file              = False,
        print_iterations                            = True,
        tol_dU                                      = 1e-3, 
        register_ref_frame                          = True,
        relax_n_iter_max                            = 30, 
        normalize_energies                          = False,
        gradient_type                               = "Sobolev",
        inner_product_H1_weight                     = 1e-3,
        mesh_degree                                 = 1,) #L2

t_stop = time.time()
print(f"duration (s) = {t_stop - t_start}")

