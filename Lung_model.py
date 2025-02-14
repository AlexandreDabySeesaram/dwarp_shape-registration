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

N_patients          = 9
Patients_Ids        = list(range(2,  N_patients + 1))

Lungs                       = ['RL','LL']


# Getting the snapshots

mapping_sphere_RL = []
mapping_sphere_LL = []
mapping_RL = []
mapping_LL = []


for patient in Patients_Ids:
    for lung in Lungs:
        if lung == "RL":
            mapping_sphere_RL.append(vtu2numpy("Results/Mapping_sphere_PA"+str(patient)+"_"+lung+"_000"))
            mapping_RL.append(vtu2numpy("Results/Mapping_PA"+str(patient)+"_"+lung+"_000"))
        else:
            mapping_sphere_LL.append(vtu2numpy("Results/Mapping_sphere_PA"+str(patient)+"_"+lung+"_000"))
            mapping_LL.append(vtu2numpy("Results/Mapping_PA"+str(patient)+"_"+lung+"_000"))

mapping_sphere_RL_array = np.hstack(mapping_sphere_RL)
mapping_sphere_LL_array = np.hstack(mapping_sphere_LL)
mapping_RL_array = np.hstack(mapping_RL)
mapping_LL_array = np.hstack(mapping_LL)

## SVD

U_sphere_RL,S_sphere_RL,V_sphere_RL = np.linalg.svd(mapping_sphere_RL_array)
U_sphere_LL,S_sphere_LL,V_sphere_LL = np.linalg.svd(mapping_sphere_LL_array)

U_RL,S_RL,V_RL = np.linalg.svd(mapping_RL_array)
U_LL,S_LL,V_LL = np.linalg.svd(mapping_LL_array)


## plot

import matplotlib.pyplot as plt

plt.figure()

plt.semilogy(S_RL, label = "Right lung")
plt.semilogy(S_LL, label = "Left lung")

plt.xlabel('Indexes of modes')
plt.ylabel('Singular values')
plt.legend()
plt.title("Lung mappings")
plt.show()
plt.close()



plt.figure()

plt.semilogy(S_sphere_RL, label = "Right lung")
plt.semilogy(S_sphere_LL, label = "Left lung")

plt.xlabel('Indexes of modes')
plt.ylabel('Singular values')
plt.legend()
plt.title("Spherical mappings")
plt.show()
plt.close()

