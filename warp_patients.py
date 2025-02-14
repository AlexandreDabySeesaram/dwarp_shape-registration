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

# Create spherical mesh


sphere_center_RL    = (120, 120, 120)                                                          
sphere_center_LL    = (320, 120, 120)                                                          
sphere_radius       = 105                                                                       
resolution          = 10                                                                                

center_RL           = dolfin.Point(sphere_center_RL[0], sphere_center_RL[1], sphere_center_RL[2])          
center_LL           = dolfin.Point(sphere_center_LL[0], sphere_center_LL[1], sphere_center_LL[2])          
radius              = sphere_radius                                                                     
domain_RL           = mshr.Sphere(center_RL, radius)
domain_LL           = mshr.Sphere(center_LL, radius)
mesh_shpere_RL      = mshr.generate_mesh(domain_RL, resolution)
mesh_shpere_LL      = mshr.generate_mesh(domain_LL, resolution)


# Import mesh

mesh_RL             = dolfin.Mesh("Meshes/mesh_RL.xml")
mesh_LL             = dolfin.Mesh("Meshes/mesh_LL.xml")

# User inputs

images_folder       = "Images/"
Results_folder      = "Results/"

N_patients          = 9
# Patients_Ids        = list(range(2,  N_patients + 1)) # Alice
Patients_Ids        = list(range(10,  40)) # Catalyn


Lungs                       = ['RL','LL']


# Pre processing

copy_pgms           = False
create_raw_vti      = False
create_bin_vti      = False
create_signed_vti   = True
from_sphere         = True

if copy_pgms:
    for patient in Patients_Ids:
        #RAW data
        source      = "/Users/daby/Documents/Code/Tracking_lungs/PA"+str(patient)+"/PGM"
        destination = "./Images/PA"+str(patient)+"/RAW"
        create_data.copy_folder(source, destination)
        #BIN data
        source      = "/Users/daby/Documents/Code/Tracking_lungs/PA"+str(patient)+"/LUNG"
        destination = "./Images/PA"+str(patient)+"/BIN"
        create_data.copy_folder(source, destination)

if create_raw_vti:
    for patient in Patients_Ids:
        create_data.PGM2vti(
            output_name         = "./Images/PA"+str(patient)+"/RAW/Pat"+str(patient)+"_RAW",
            Raw_PGM_base_name   = "./Images/PA"+str(patient)+"/RAW/Pat"+str(patient)+"_inspi1",
                )


if create_raw_vti:
    for patient in Patients_Ids:
        create_data.PGM2vti(
            output_name         = "./Images/PA"+str(patient)+"/RAW/Pat"+str(patient)+"_RAW",
            Raw_PGM_base_name   = "./Images/PA"+str(patient)+"/RAW/Pat"+str(patient)+"_inspi1",
                )

if create_bin_vti:
    for patient in Patients_Ids:
        create_data.PGM2vti(
            output_name         = "./Images/PA"+str(patient)+"/BIN/Pat"+str(patient)+"_BIN",
            Raw_PGM_base_name   = "./Images/PA"+str(patient)+"/RAW/Pat"+str(patient)+"_inspi1",
            Bin_PGM_base_name   = "./Images/PA"+str(patient)+"/BIN/Pat"+str(patient)+"_inspi1",
                )

if create_signed_vti:

    for patient in Patients_Ids:
        for lung in Lungs:
            if lung == "RL":
                create_data.sign_masking_binary(
                    input_name          = "./Images/PA"+str(patient)+"/BIN/Pat"+str(patient)+"_BIN",
                    suffix              = "signed_RL_0"                  , 
                    scalar2zero         = 200                           ,
                    scalar_background   = 0                             ,                               # Initial background pixel intensity
                    scalar_foreground   = 100                           ,                               # Initial foreground pixel intensity
                    target_value_bg     = 1                             ,                               # target background pixel intensity
                    target_value_fg     = -1                            ,                               # target foreground pixel intensity
                    target_type         = "signed_char"                 ,                               # unsigned_char, signed_char, float
                    )
            elif lung == "LL":
                create_data.sign_masking_binary(
                    input_name          = "./Images/PA"+str(patient)+"/BIN/Pat"+str(patient)+"_BIN",
                    suffix              = "signed_LL_0"                  , 
                    scalar2zero         = 100                           ,
                    scalar_background   = 0                             ,                               # Initial background pixel intensity
                    scalar_foreground   = 200                           ,                               # Initial foreground pixel intensity
                    target_value_bg     = 1                             ,                               # target background pixel intensity
                    target_value_fg     = -1                            ,                               # target foreground pixel intensity
                    target_type         = "signed_char"                 ,                               # unsigned_char, signed_char, float
                    )




## Tracking

result_folder               = "./Results"
reduced_kinematics_model    = "translation+scaling"

for patient in Patients_Ids:
    images_folder   = "./Images/PA"+str(patient)+"/BIN"
    for lung in Lungs:
        if lung == "RL":
            image_name      = "Pat"+str(patient)+"_BIN"+"_signed_RL"
            # mesh            = mesh_RL
            if from_sphere:
                mesh            = mesh_shpere_RL
            else:
                mesh            = mesh_RL
        elif lung == "LL":
            image_name      = "Pat"+str(patient)+"_BIN"+"_signed_LL"
            # mesh            = mesh_LL
            if from_sphere:
                mesh            = mesh_shpere_LL
            else:
                mesh            = mesh_LL

        if from_sphere:
            filebasename    = "Mapping_sphere_PA"+str(patient)+"_"+lung
        else:
            filebasename    = "Mapping_PA"+str(patient)+"_"+lung

        processing.reduced_kinematics(        
            result_folder                          = result_folder              ,
            filebasename                           = filebasename               ,
            image_folder                           = images_folder              ,
            image_name                             = image_name                 ,
            mesh                                   = mesh                       ,
            reduced_kinematics_model               = reduced_kinematics_model   ,
            tol_dU                                 = 1e-3                       ,
            n_iter_max                             = 100                        ,
            images_quadrature                      = 6                          ,
            )

        processing.full_kinematics(        
            result_folder                          = result_folder              ,
            filebasename                           = filebasename               ,
            image_folder                           = images_folder              ,
            image_name                             = image_name                 ,
            mesh                                   = mesh                       ,
            tol_dU                                 = 1e-3                       ,
            n_iter_max                             = 100                        ,
            images_quadrature                      = 6                          ,
            )