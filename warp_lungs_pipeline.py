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
resolution      = 35                                                                        # Resolution of the mesh

# Create a 3D spherical domain
center          = dolfin.Point(sphere_center[0], sphere_center[1], sphere_center[2])        # Center of the disc
radius          = 1.4*sphere_radius                                                         # Radius of the disc
domain          = mshr.Sphere(center, radius)
mesh_shpere_RL      = mshr.generate_mesh(domain, resolution)
mesh_shpere_RL.num_vertices()


# sphere_center_RL    = (120, 120, 120)                                                          
# sphere_center_LL    = (310, 120, 120)                                                          
# sphere_radius       = 105                                                                       
# resolution          = 10                                                                                

# center_RL           = dolfin.Point(sphere_center_RL[0], sphere_center_RL[1], sphere_center_RL[2])          
# center_LL           = dolfin.Point(sphere_center_LL[0], sphere_center_LL[1], sphere_center_LL[2])          
# radius              = sphere_radius                                                                     
# domain_RL           = mshr.Sphere(center_RL, radius)
# domain_LL           = mshr.Sphere(center_LL, radius)
# mesh_shpere_RL      = mshr.generate_mesh(domain_RL, resolution)
# mesh_shpere_LL      = mshr.generate_mesh(domain_LL, resolution)
# print(f"LL number of nodes: {mesh_shpere_LL.num_vertices()}")
# print(f"RL number of nodes: {mesh_shpere_RL.num_vertices()}")

# # Import mesh

# mesh_RL             = dolfin.Mesh("Meshes/mesh_RL.xml")
# mesh_LL             = dolfin.Mesh("Meshes/mesh_LL.xml")

#%% load previous lung mesh


mesh_lung=dolfin.Mesh("Meshes/mesh_RL.xml")


#%% Load image

mesh_name                   = "3D_lung_PA5"
image_basename              = "Pat2_BIN"
image_suffix                = "signed_RL"
result_folder               = "Results" 
filebasename                = "test_mapping_PA2_L2"
image_name                  = image_basename+"_"+image_suffix
image_folder                = "Images"
image_path                  = image_folder+"/"+image_name
reduced_kinematics_model    = "translation+scaling"


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
        working_basename                            = filebasename+"_reduced",
        images_folder                               = image_folder,
        images_basename                             = image_name,
        images_ext                                  = "vti",
        mesh                                        = mesh_shpere_RL,
        kinematics_type                             = "reduced",
        reduced_kinematics_model                    = reduced_kinematics_model,
        images_quadrature                           = 6,
        n_iter_max                                  = 100,
        nonlinearsolver                             = "gradient_descent",
        min_gradient_step                           = 1e-6,
        gradient_step                               = 1e-2,  
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
        inner_product_H1_weight                     = 1e-3) #L2

t_stop = time.time()
print(f"duration (s) = {t_stop - t_start}")


dwarp.warp(
        warping_type                                = "registration"                                            , # registration
        working_folder                              = result_folder,
        images_char_func                            = False,
        working_basename                            = filebasename,
        images_folder                               = image_folder,
        images_basename                             = image_name,
        images_ext                                  = "vti",
        mesh                                        = mesh_shpere_RL,
        kinematics_type                             = "full",
        reduced_kinematics_model                    = reduced_kinematics_model,
        nonlinearsolver                             = "gradient_descent",
        images_quadrature                           = 6,
        n_iter_max                                  = 100,
        # relax_type                                  = "gss",
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
        initialize_U_from_file                      = True,
        initialize_U_folder                         = result_folder, 
        initialize_U_basename                       = filebasename+"_reduced",
        initialize_U_ext                            = "vtu",
        initialize_U_array_name                     = "displacement",
        initialize_U_method                         = "interpolation") # dofs_transfer, interpolation, projection) #L2