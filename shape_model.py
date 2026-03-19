import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import pyvista as pv
import matplotlib.pyplot as plt


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

# Load the VTU file
vtu_name = "Results/Mapping_sphere_PA"+str(12)+"_RL_000"  # Replace with your file name
_ , output = vtu2numpy(vtu_name, return_output = True)


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


#%% Stack the 3 space components

reduced_mapping_sphere_RL_stacked = [mapping.flatten()[:,None] for mapping in reduced_mapping_sphere_RL]
reduced_mapping_sphere_LL_stacked = [mapping.flatten()[:,None] for mapping in reduced_mapping_sphere_LL]

mapping_sphere_RL_stacked = [mapping.flatten()[:,None] for mapping in mapping_sphere_RL]
mapping_sphere_LL_stacked = [mapping.flatten()[:,None] for mapping in mapping_sphere_LL]

if shperical_mapping:
    mapping_sphere_RL_array = np.hstack(mapping_sphere_RL_stacked)
    mapping_sphere_LL_array = np.hstack(mapping_sphere_LL_stacked)
    reduced_mapping_sphere_RL_array = np.hstack(reduced_mapping_sphere_RL_stacked)
    reduced_mapping_sphere_LL_array = np.hstack(reduced_mapping_sphere_LL_stacked)


    U_soft_sphere_RL,S_soft_sphere_RL,V_soft_sphere_RL = np.linalg.svd(mapping_sphere_RL_array - reduced_mapping_sphere_RL_array)
    U_soft_sphere_LL,S_soft_sphere_LL,V_soft_sphere_LL = np.linalg.svd(mapping_sphere_LL_array - reduced_mapping_sphere_LL_array)

if lung_mapping:
    mapping_RL_array        = np.hstack(mapping_RL)
    mapping_LL_array        = np.hstack(mapping_LL)
    reduced_mapping_RL_array        = np.hstack(reduced_mapping_RL)
    reduced_mapping_LL_array        = np.hstack(reduced_mapping_LL)




plt.figure()

plt.semilogy(S_soft_sphere_RL[:-4], label = "Right lung")
plt.semilogy(S_soft_sphere_LL[:-4], label = "Left lung")

plt.xlabel('Indexes of modes')
plt.ylabel('Singular values')
plt.legend()
plt.title("Spherical mappings full mapping")
plt.show()
plt.close()




#%% Fixed point

soft_mapping_flat = mapping_sphere_RL_array - reduced_mapping_sphere_RL_array

np.save("soft_mappings_flat.npy", soft_mapping_flat)


count = 0 
eta = 1.2
# while eta >=1e-5 and count <= 500:
#     count += 1
#     soft_mapping_flat_old = soft_mapping_flat.copy()
#     soft_mapping_array = soft_mapping_flat.reshape((-1,3,39))
#     centers_mass = soft_mapping_array.sum(axis=0)/soft_mapping_array.shape[0]
#     soft_mapping_array -= centers_mass[None, :, :]                              # centered mappings

#     soft_mapping_array_average = soft_mapping_array.sum(axis=2)/soft_mapping_array.shape[2]
#     soft_mapping_array -= soft_mapping_array_average[ :, :, None]                              # remove average
#     soft_mapping_flat = soft_mapping_array.reshape((-1,39))
#     eta = np.mean(np.linalg.norm(soft_mapping_flat - soft_mapping_flat_old, axis=0)/np.linalg.norm(soft_mapping_flat, axis=0))
#     print(f"* eta = {eta}, count = {count}")






U_soft_sphere,S_soft_sphere,V_soft_sphere = np.linalg.svd(soft_mapping_flat)

mean_average = (soft_mapping_flat.reshape((-1,3,39))[:,:3,0] + soft_mapping_flat.reshape((-1,3,39))[:,:3,1])/2

#%% Plot
# Convert the VTK object to a PyVista mesh
mesh = pv.wrap(output)

saved_initial_mesh = mesh.points.copy()

# Plot the modified mesh interactively
plotter = pv.Plotter()
plotter.add_mesh(mesh, show_scalar_bar=True)
plotter.show()

# Manually wrap the mesh by the displacement field
avrg_flat = np.load("Average_mapping.npy")
avrg = avrg_flat.reshape((-1,3))
mesh.points += avrg

# Update the plot
plotter.render()

