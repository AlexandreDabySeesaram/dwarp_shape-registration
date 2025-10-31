import dolfin
import dolfin_warp as dwarp
import os
import glob



def compute_barycenter(
                    mesh,
                    lung                        = "RL", 
                    basename                    = "Barycenter_fine_mapping", 
                    mappings_basename           = "Mapping_fine_sphere", 
                    regul_model                 = "ogdenciarletgeymonatneohookean",
                    relax_type                  = "aitken", 
                    initialize_U_from_file      = True                                ,
                    initialize_U_folder         = "Initialisation"                    ,
                    initialize_U_basename       = "Mapping_coarse_sphere_RL"          ,
                    initialize_U_ext            = "vtu"                               ,
                    initialize_U_array_name     = "displacement"                      ,
                    initialize_U_method         = "dofs_transfer"                     ,   ):
    dwarp.warp(
            working_folder                              = "Results/barycenter", 
            working_basename                            = basename+"_"+lung,
            images_folder                               = "None",
            images_basename                             = "None",
            mappings_folder                             = "Mappings",
            # mappings_basename                           = "Mapping_fine_sphere_interpolated_RL",
            mappings_basename                           = mappings_basename+"_"+lung,
            mesh                                        = mesh,
            warping_type                                = "barycenter"                                            , # registration
            kinematics_type                             = "full",
            regul_model                                 = regul_model, # ogdenciarletgeymonatneohookean
            # nonlinearsolver                             = "newton",
            nonlinearsolver                             = "gradient_descent",
            images_quadrature                           = 6,
            n_iter_max                                  = 10000,
            relax_type                                  = relax_type,
            min_gradient_step                           = 1e-6,
            gradient_step                               = 1,
            write_VTU_files                             = True,
            write_VTU_files_with_preserved_connectivity = True,
            print_iterations                            = True,
            tol_dU                                      = 1e-6, 
            relax_n_iter_max                            = 30, 
            normalize_energies                          = False,
            gradient_type                               = "L2",
            continue_after_fail                         = True,
            inner_product_H1_weight                     = 1e-3, 
            initialize_U_from_file                      = initialize_U_from_file                                ,
            initialize_U_folder                         = initialize_U_folder                                ,
            initialize_U_basename                       = initialize_U_basename                                ,
            initialize_U_ext                            = initialize_U_ext                              ,
            initialize_U_array_name                     = initialize_U_array_name                  ,
            initialize_U_method                         = initialize_U_method                     , # dofs_transfer, interpolation, projection

            )




mesh_coarse_RL = dolfin.Mesh("Meshes/Coarse_sphere_RL.xml")
mesh_coarse_LL = dolfin.Mesh("Meshes/Coarse_sphere_LL.xml")
mesh_fine_RL   = dolfin.Mesh("Meshes/Fine_sphere_RL.xml")
mesh_fine_LL   = dolfin.Mesh("Meshes/Fine_sphere_LL.xml")



model               = "ogdenciarletgeymonatneohookean"          # ogdenciarletgeymonatneohookean, hooke
lung                = "RL"
coarsness           = "fine"
basename            = "init_2_Barycenter_"+ coarsness+"_"+model
# mappings_basename   = "Mapping_"+coarsness+"_sphere"
mappings_basename   = "Mapping_"+coarsness+"_sphere"

match lung:
    case "RL": 
        if coarsness == "coarse":
            mesh = mesh_coarse_RL
        else:
            mesh = mesh_fine_RL
    case "LL": 
        if coarsness == "coarse":
            mesh = mesh_coarse_LL
        else:
            mesh = mesh_fine_LL

for vtu_filename in glob.glob("Results/barycenter/"+basename+"_"+lung+"-frame=None"+"_[0-9]*.vtu"):
    os.remove(vtu_filename)

 

compute_barycenter(
                    mesh,
                    lung                        = lung, 
                    basename                    = basename, 
                    mappings_basename           = mappings_basename, 
                    regul_model                 = model,
                    relax_type                  = "backtracking", 
                    initialize_U_from_file      = True                                ,
                    initialize_U_folder         = "Initialisation"                    ,
                    initialize_U_basename       = "Mapping_"+coarsness+"_sphere_RL"   ,
                    initialize_U_ext            = "vtu"                               ,
                    initialize_U_array_name     = "displacement"                      ,
                    initialize_U_method         = "dofs_transfer"                     , 

                        ) 