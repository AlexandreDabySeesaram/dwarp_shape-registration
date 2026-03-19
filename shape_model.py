import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import pyvista as pv
import matplotlib.pyplot as plt




job = "configurations/barycenter_arm_config"
try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # Older versions

with open(f"{job}.toml", "rb") as f:
    config = tomllib.load(f)



def update_plot(mesh, mapping, initial_nodes):

    mesh.points = initial_nodes.copy()
    mesh.points += mapping

    plotter.render()

def vtu2numpy(vtu_name, return_output = False):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_name + ".vtu")
    print(f"Name : {vtu_name + '.vtu'}")
    reader.Update()
    output = reader.GetOutput()
    displacement_vtk_array = output.GetPointData().GetArray('displacement')
    displacement_array = vtk_to_numpy(displacement_vtk_array)
    if return_output:
        return displacement_array, output
    else:
        return displacement_array


# user input

N_patients                      = 40
Patients_Ids                    = list(range(2,  N_patients + 1))


Lungs                           = ['RL']


# Getting the snapshots

mapping_RL           = []
mapping_LL           = []
reduced_mapping_RL   = [] 
reduced_mapping_LL   = []
reduced_mapping_LL          = [] 


for patient in Patients_Ids:
    for lung in Lungs:
        if lung == "RL": 
          mapping_RL.append(vtu2numpy("Results/Mapping_PA"+str(patient)+"_"+lung+"_000"))
          reduced_mapping_RL.append(vtu2numpy("Results/Mapping_PA"+str(patient)+"_"+lung+"_reduced_000"))
        
        else:
          mapping_LL.append(vtu2numpy("Results/Mapping_PA"+str(patient)+"_"+lung+"_000"))
          reduced_mapping_LL.append(vtu2numpy("Results/Mapping_PA"+str(patient)+"_"+lung+"_reduced_000"))
            
#%% Stack the 3 space components

reduced_mapping_RL_stacked = [mapping.flatten()[:,None] for mapping in reduced_mapping_RL]
reduced_mapping_LL_stacked = [mapping.flatten()[:,None] for mapping in reduced_mapping_LL]

mapping_RL_stacked = [mapping.flatten()[:,None] for mapping in mapping_RL]
mapping_LL_stacked = [mapping.flatten()[:,None] for mapping in mapping_LL]

mapping_RL_array = np.hstack(mapping_RL_stacked)
mapping_LL_array = np.hstack(mapping_LL_stacked)
reduced_mapping_RL_array = np.hstack(reduced_mapping_RL_stacked)
reduced_mapping_LL_array = np.hstack(reduced_mapping_LL_stacked)


U_soft_RL,S_soft_RL,V_soft_RL = np.linalg.svd(mapping_RL_array - reduced_mapping_RL_array)
    

plt.figure()

plt.semilogy(S_soft_RL[:-4], label = "Right lung")
plt.semilogy(S_soft_LL[:-4], label = "Left lung")

plt.xlabel('Indexes of modes')
plt.ylabel('Singular values')
plt.legend()
plt.title("Spherical mappings full mapping")
plt.show()
plt.close()




