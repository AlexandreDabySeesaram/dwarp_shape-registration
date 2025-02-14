import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np

def vtu2numpy(vtu_name):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_name+".vtu")

    print(f"Name : {vtu_name+'.vtu'}")

    reader.Update()

    output = reader.GetOutput()


    displacement_vtk_array = output.GetPointData().GetArray('displacement')

    displacement_array = vtk_to_numpy(displacement_vtk_array)
    displacement_array_flat = displacement_array.flatten().reshape(-1,1)
    return displacement_array_flat


# user input

N_patients                      = 9
# Patients_Ids                    = list(range(2,  N_patients + 1))
Patients_Ids                    = list(range(2,  41))

Lungs                           = ['RL','LL']

shperical_mapping               = True
lung_mapping                    = False

# Getting the snapshots

if shperical_mapping:
    mapping_sphere_RL           = []
    mapping_sphere_LL           = []
    reduced_mapping_sphere_RL   = []
    reduced_mapping_sphere_LL   = []
if lung_mapping:  
    mapping_RL                  = []
    mapping_LL                  = []
    reduced_mapping_RL          = []
    reduced_mapping_LL          = []


for patient in Patients_Ids:
    for lung in Lungs:
        if lung == "RL":
            if shperical_mapping:
                mapping_sphere_RL.append(vtu2numpy("Results/Mapping_sphere_PA"+str(patient)+"_"+lung+"_000"))
                reduced_mapping_sphere_RL.append(vtu2numpy("Results/Mapping_sphere_PA"+str(patient)+"_"+lung+"_reduced_000"))
            if lung_mapping:
                mapping_RL.append(vtu2numpy("Results/Mapping_PA"+str(patient)+"_"+lung+"_000"))
                reduced_mapping_RL.append(vtu2numpy("Results/Mapping_PA"+str(patient)+"_"+lung+"_reduced_000"))
        else:
            if shperical_mapping:
                mapping_sphere_LL.append(vtu2numpy("Results/Mapping_sphere_PA"+str(patient)+"_"+lung+"_000"))
                reduced_mapping_sphere_LL.append(vtu2numpy("Results/Mapping_sphere_PA"+str(patient)+"_"+lung+"_reduced_000"))
            if lung_mapping:
                mapping_LL.append(vtu2numpy("Results/Mapping_PA"+str(patient)+"_"+lung+"_000"))
                reduced_mapping_LL.append(vtu2numpy("Results/Mapping_PA"+str(patient)+"_"+lung+"_reduced_000"))

if shperical_mapping:
    mapping_sphere_RL_array = np.hstack(mapping_sphere_RL)
    mapping_sphere_LL_array = np.hstack(mapping_sphere_LL)
    reduced_mapping_sphere_RL_array = np.hstack(reduced_mapping_sphere_RL)
    reduced_mapping_sphere_LL_array = np.hstack(reduced_mapping_sphere_LL)
if lung_mapping:
    mapping_RL_array        = np.hstack(mapping_RL)
    mapping_LL_array        = np.hstack(mapping_LL)
    reduced_mapping_RL_array        = np.hstack(reduced_mapping_RL)
    reduced_mapping_LL_array        = np.hstack(reduced_mapping_LL)

## SVD
if shperical_mapping:
    U_sphere_RL,S_sphere_RL,V_sphere_RL = np.linalg.svd(mapping_sphere_RL_array)
    U_sphere_LL,S_sphere_LL,V_sphere_LL = np.linalg.svd(mapping_sphere_LL_array)
    U_soft_sphere_RL,S_soft_sphere_RL,V_soft_sphere_RL = np.linalg.svd(mapping_sphere_RL_array - reduced_mapping_sphere_RL_array)
    U_soft_sphere_LL,S_soft_sphere_LL,V_soft_sphere_LL = np.linalg.svd(mapping_sphere_LL_array - reduced_mapping_sphere_LL_array)

if lung_mapping:
    U_RL,S_RL,V_RL = np.linalg.svd(mapping_RL_array)
    U_LL,S_LL,V_LL = np.linalg.svd(mapping_LL_array)
    U_soft_RL,S_soft_RL,V_soft_RL = np.linalg.svd(mapping_RL_array - reduced_mapping_RL_array)
    U_soft_LL,S_soft_LL,V_soft_LL = np.linalg.svd(mapping_LL_array - reduced_mapping_LL_array)


## plot

import matplotlib.pyplot as plt


if lung_mapping:

    plt.figure()

    plt.semilogy(S_RL, label = "Right lung")
    plt.semilogy(S_LL, label = "Left lung")

    plt.xlabel('Indexes of modes')
    plt.ylabel('Singular values')
    plt.legend()
    plt.title("Lung mappings")
    plt.show()
    plt.close()


if shperical_mapping:

    plt.figure()

    plt.semilogy(S_sphere_RL, label = "Right lung")
    plt.semilogy(S_sphere_LL, label = "Left lung")

    plt.xlabel('Indexes of modes')
    plt.ylabel('Singular values')
    plt.legend()
    plt.title("Spherical mappings full mapping")
    plt.show()
    plt.close()

    plt.figure()

    plt.semilogy(S_soft_sphere_RL[:-4], label = "Right lung")
    plt.semilogy(S_soft_sphere_LL[:-4], label = "Left lung")

    plt.xlabel('Indexes of modes')
    plt.ylabel('Singular values')
    plt.legend()
    plt.title("Spherical mappings without reduced-kinematics")
    plt.show()
    plt.close()

import pandas as pd

x = np.linspace(0,S_soft_sphere_RL[:-4].shape[0],S_soft_sphere_RL[:-4].shape[0])
y_RL = S_soft_sphere_RL[:-4]
y_LL = S_soft_sphere_LL[:-4]
y_RL_full = S_sphere_RL[:-4]
y_LL_full = S_sphere_LL[:-4]

df = pd.DataFrame(np.stack((x,y_RL,y_LL,y_RL_full,y_LL_full)).T, columns=['N','SV_RL','SV_LL','SV_full_RL','SV_full_LL'])
df.to_csv('Results/Figure_Singular_values.csv', index=False)

