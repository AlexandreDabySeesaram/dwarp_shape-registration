import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import dolfin
import dolfin_warp as dwarp
import mshr
import dolfin_mech as dmech
import create_data
import glob
import processing

job = "configurations/registration_config"
try:
    import tomllib
    config = tomllib.load(f)
except: 
    import tomli
    with open(job+".toml", "rb") as f:
        config = tomli.load(f)


mesh_dict = {}

match config["tracking"]["mesh_initialisation"]:
    case "sphere":
        match config["tracking"]["fineness"]:
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

        # print(f"RL mesh generated - number of nodes: {mesh_shpere_RL.num_vertices()}")
        # print(f"LL mesh generated - number of nodes: {mesh_shpere_LL.num_vertices()}")

    case "lung":
        if "RL" in config["tracking"]["lungs"]:
            mesh_dict["RL"] = dolfin.Mesh(config["names"]["mesh_folder"]+"/mesh_RL.xml")
            print("* RL mesh loaded from lung mesh")
        if "LL" in config["tracking"]["lungs"]:
            mesh_dict["LL"] = dolfin.Mesh(config["names"]["mesh_folder"]+"/mesh_LL.xml")
            print("* LL mesh loaded from lung mesh")

    case "morphed_sphere":
        if "RL" in config["tracking"]["lungs"]:
            mesh_dict["RL"] = dolfin.Mesh(config["names"]["mesh_folder"]+"/"+config["tracking"]["fineness"]+"_morphed_sphere_RL.xml")
            print("* RL mesh loaded from morphed_sphere")
        if "LL" in config["tracking"]["lungs"]:
            mesh_dict["LL"] = dolfin.Mesh(config["names"]["mesh_folder"]+"/"+config["tracking"]["fineness"]+"_morphed_sphere_LL.xml")
            print("* LL mesh loaded from morphed_sphere")


Patients_Ids        = list(range(config["tracking"]["patients"][0],  config["tracking"]["patients"][1])) 



## Tracking

result_folder               = config["names"]["result_folder"]
reduced_kinematics_model    = config["tracking"]["reduced_kinematics_model"]

for patient in Patients_Ids:
    images_folder   = config["names"]["image_folder"]+"/PA"+str(patient)+"/BIN"
    for lung in config["tracking"]["lungs"]:
        image_name      = "Pat"+str(patient)+"_BIN"+"_signed_"+lung
        mesh            = mesh_dict[lung]
        filebasename    = config["names"]["filebasename"]+"_"+config["tracking"]["fineness"]+"_mesh_"+config["tracking"]["mesh_initialisation"]+"_PA"+str(patient)+"_"+lung

        processing.reduced_kinematics(        
            result_folder                          = result_folder                                  ,
            filebasename                           = filebasename                                   ,
            image_folder                           = images_folder                                  ,
            image_name                             = image_name                                     ,
            mesh                                   = mesh                                           ,
            reduced_kinematics_model               = config["tracking"]["reduced_kinematics_model"] ,
            tol_dU                                 = config["tracking"]["tol_dU_RK"]                ,
            n_iter_max                             = config["tracking"]["n_iter_max_RK"]            ,
            images_quadrature                      = config["tracking"]["images_quadrature_RK"]     ,
            )

        processing.full_kinematics(        
            result_folder                          = result_folder                                  ,
            filebasename                           = filebasename                                   ,
            image_folder                           = images_folder                                  ,
            image_name                             = image_name                                     ,
            mesh                                   = mesh                                           ,
            tol_dU                                 = config["tracking"]["tol_dU_RK"]                ,
            n_iter_max                             = config["tracking"]["n_iter_max_RK"]            ,
            images_quadrature                      = config["tracking"]["images_quadrature_RK"]     ,
            )