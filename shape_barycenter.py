import dolfin
import dolfin_warp as dwarp
import os
import glob

#job = "configurations/registration_config"
job = "configurations/barycenter_config"
try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # Older versions

with open(f"{job}.toml", "rb") as f:
    config = tomllib.load(f)


def compute_barycenter(
        mesh,
        lung                        = "RL"                                     , 
        basename                    = "Barycenter_fine_mapping"                , 
        mappings_basename           = "Mapping_fine_sphere"                    , 
        regul_model                 = "ogdenciarletgeymonatneohookean"         ,
        relax_type                  = "aitken"                                 , 
        initialize_U_from_file      = True                                     ,
        initialize_U_folder         = "Initialisation"                         ,
        initialize_U_basename       = "Mapping_coarse_sphere_RL"               ,
        initialize_U_ext            = "vtu"                                    ,
        initialize_U_array_name     = "displacement"                           ,
        initialize_U_method         = "dofs_transfer"                          ,):
    dwarp.warp(
            working_folder                                 = "Results/barycenter"          , 
            working_basename                               = basename+"_"+lung             ,
            images_folder                                  = "None"                        ,
            images_basename                                = "None"                        ,
            mappings_folder                                = "Mappings"                    ,
            mappings_basename                              = mappings_basename+"_"+lung    ,
            mesh                                           = mesh                          ,
            warping_type                                   = "barycenter"                  , # registration
            kinematics_type                                = "full"                        ,
            regul_model                                    = regul_model                   , # ogdenciarletgeymonatneohookean
            nonlinearsolver                                = "gradient_descent"            ,
            images_quadrature                              = 6                             ,
            n_iter_max                                     = 10000                         ,
            relax_type                                     = relax_type                    ,
            min_gradient_step                              = 1e-6                          ,
            gradient_step                                  = 1                             ,    
            write_VTU_files                                = True                          ,
            write_VTU_files_with_preserved_connectivity    = True                          ,
            print_iterations                               = True                          ,
            tol_dU                                         = 1e-6                          , 
            relax_n_iter_max                               = 30                            , 
            normalize_energies                             = False                         ,
            gradient_type                                  = "L2"                          ,
            continue_after_fail                            = True                          ,
            inner_product_H1_weight                        = 1e-3                          , 
            initialize_U_from_file                         = initialize_U_from_file        ,
            initialize_U_folder                            = initialize_U_folder           ,
            initialize_U_basename                          = initialize_U_basename         ,
            initialize_U_ext                               = initialize_U_ext              ,
            initialize_U_array_name                        = initialize_U_array_name       ,
            initialize_U_method                            = initialize_U_method           , # dofs_transfer, interpolation, projection

            )


 
mesh_dict = {}

match config["mappings"]["mesh_initialisation"]:
    case "sphere":
        match config["mappings"]["fineness"]:
            case "fine":
                resolution = 35
            case "Coarse":
                resolution = 5

        sphere_center_RL    = (120, 120, 120)                                                          
        sphere_center_LL    = (310, 120, 120)                                                          
        sphere_radius       = 105                                                                       
        center_RL           = dolfin.Point(sphere_center_RL[0], sphere_center_RL[1], sphere_center_RL[2])          
        center_LL           = dolfin.Point(sphere_center_LL[0], sphere_center_LL[1], sphere_center_LL[2])          
        radius              = sphere_radius                                                                     
        domain_RL           = mshr.Sphere(center_RL, radius)
        domain_LL           = mshr.Sphere(center_LL, radius)
        mesh_dict["RL"]     = mshr.generate_mesh(domain_RL, resolution)
        mesh_dict["LL"]     = mshr.generate_mesh(domain_LL, resolution)


    case "lung":
        if "RL" in config["range"]["lungs"]:
            mesh_dict["RL"] = dolfin.Mesh(config["names"]["mesh_folder"]+"/mesh_RL.xml")
            print("* RL mesh loaded from lung mesh")
        if "LL" in config["range"]["lungs"]:
            mesh_dict["LL"] = dolfin.Mesh(config["names"]["mesh_folder"]+"/mesh_LL.xml")
            print("* LL mesh loaded from lung mesh")

    case "morphed_sphere":
        if "RL" in config["range"]["lungs"]:
            mesh_dict["RL"] = dolfin.Mesh(config["names"]["mesh_folder"]+"/"+config["mappings"]["fineness"]+"_morphed_sphere_RL.xml")
            print("* RL mesh loaded from morphed_sphere")
        if "LL" in config["range"]["lungs"]:
            mesh_dict["LL"] = dolfin.Mesh(config["names"]["mesh_folder"]+"/"+config["mappings"]["fineness"]+"_morphed_sphere_LL.xml")
            print("* LL mesh loaded from morphed_sphere")


Patients_Ids        = list(range(config["range"]["patients"][0],  config["range"]["patients"][1])) 

model               = config["barycenter"]["model"]
lungs               = config["range"]["lungs"]
mappings_basename   = config["names"]["mappings_basename"]+config["mappings"]["fineness"]+config["mappings"]["mesh_initialisation"]
basename            = config["barycenter"]["basename"]


for lung in lungs:

    for vtu_filename in glob.glob("Results/barycenter/"+basename+"_"+lung+"-frame=None"+"_[0-9]*.vtu"):
        os.remove(vtu_filename)

    compute_barycenter(
            mesh                        = mesh_dict[lung],
            lung                        = lung, 
            basename                    = basename+lung, 
            mappings_basename           = mappings_basename, 
            regul_model                 = model,
            relax_type                  = "backtracking", 
            initialize_U_from_file      = False                               ,
            # initialize_U_folder         = "Initialisation"                    ,
            # initialize_U_basename       = "Mapping_"+coarsness+"_sphere_RL"   ,
            # initialize_U_ext            = "vtu"                               ,
            # initialize_U_array_name     = "displacement"                      ,
            # initialize_U_method         = "dofs_transfer"                     , 
            ) 
