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




def plot_latent_space_3d(S, Vh):
    """
    S: Singular values (S_soft_RL_centered)
    Vh: Right singular vectors (V_soft_RL_centered)
    """
    # 1. Calculate the Latent Coordinates (Scores)
    # We take the first 3 components. 
    # Coordinates = S_i * Vh_i
    # This gives us the position of each 'sample' in the 3D reduced space.
    pc1 = S[0] * Vh[0, :]
    pc2 = S[1] * Vh[1, :]
    pc3 = S[2] * Vh[2, :]

    # 2. Create the 3D Plot
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot
    scatter = ax.scatter(pc1, pc2, pc3, c=pc1, cmap='viridis', s=50, alpha=0.8)

    # Labeling
    ax.set_title('3D Latent Space (Top 3 Principal Components)', fontsize=15)
    ax.set_xlabel(f'PC1 (Var: {S[0]**2 / np.sum(S**2):.1%})')
    ax.set_ylabel(f'PC2 (Var: {S[1]**2 / np.sum(S**2):.1%})')
    ax.set_zlabel(f'PC3 (Var: {S[2]**2 / np.sum(S**2):.1%})')

    # Add a colorbar to show the gradient of the first component
    fig.colorbar(scatter, ax=ax, label='PC1 Magnitude', shrink=0.5)

    plt.tight_layout()
    plt.savefig('latent_space_3d.png')
    plt.show()













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

mapping_base_name = "Mappings/Mapping_arm_Fine_mesh_morphed_sphere"

for patient in Patients_Ids:
    for lung in Lungs:
        if lung == "RL": 
            try:    
                mapping_RL.append(vtu2numpy(mapping_base_name+"_"+lung+"_"+str(patient)))
                reduced_mapping_RL.append(vtu2numpy(mapping_base_name+"_reduced"+"_"+lung+"_"+str(patient)))
            except:
                print("Mapping not in database")
        else:
            try:    
                mapping_LL.append(vtu2numpy(mapping_base_name+"_"+lung)+str(patient))
                reduced_mapping_LL.append(vtu2numpy(mapping_base_name+"_reduced"+"_"+lung)+str(patient))
            except:
                print("Mapping not in database")
            
#%% Stack the 3 space components

barycenter_RL = vtu2numpy("Results/barycenter/barycenter_arm_Fine_mesh_morphed_sphere_ogdenciarletgeymonatneohookean_RL-frame=None_5130")
barycenter_RL_stacked = barycenter_RL.flatten()[:,None]

reduced_mapping_RL_stacked = [mapping.flatten()[:,None] for mapping in reduced_mapping_RL]
# reduced_mapping_LL_stacked = [mapping.flatten()[:,None] for mapping in reduced_mapping_LL]

mapping_RL_stacked = [mapping.flatten()[:,None] for mapping in mapping_RL]
# mapping_LL_stacked = [mapping.flatten()[:,None] for mapping in mapping_LL]

mapping_RL_array = np.hstack(mapping_RL_stacked)
# mapping_LL_array = np.hstack(mapping_LL_stacked)
reduced_mapping_RL_array = np.hstack(reduced_mapping_RL_stacked)
# reduced_mapping_LL_array = np.hstack(reduced_mapping_LL_stacked)


U_soft_RL,S_soft_RL,V_soft_RL = np.linalg.svd(mapping_RL_array - reduced_mapping_RL_array)
    

plt.figure()

plt.semilogy(S_soft_RL[:-4], label = "Right lung")
# plt.semilogy(S_soft_LL[:-4], label = "Left lung")

plt.xlabel('Indexes of modes')
plt.ylabel('Singular values')
plt.legend()
plt.title("Spherical mappings full mapping")
plt.show()
plt.close()


U_soft_RL_centered,S_soft_RL_centered,V_soft_RL_centered = np.linalg.svd(mapping_RL_array - reduced_mapping_RL_array - barycenter_RL_stacked)
    

plt.figure()

plt.semilogy(S_soft_RL_centered[:-4], label = "Right lung")
# plt.semilogy(S_soft_LL[:-4], label = "Left lung")

plt.xlabel('Indexes of modes')
plt.ylabel('Singular values')
plt.legend()
plt.title("Spherical mappings full mapping")
plt.show()
plt.close()





# Assuming you already have:
# U_soft_RL_centered, S_soft_RL_centered, V_soft_RL_centered
# barycenter_RL_stacked

def reconstruct_shape(U, S, barycenter, slider_vals):
    """
    U: Left singular vectors from SVD
    S: Singular values
    barycenter: The mean shape used for centering
    slider_vals: A list/array of 3 multipliers (e.g., [1.0, -0.5, 2.0])
    """
    # 1. Start with the mean shape (barycenter)
    # We copy it so we don't overwrite the original
    new_shape = np.copy(barycenter)
    
    # 2. Add the contribution of the first 3 components
    # Each component is U[:, i] * S[i] * slider_multiplier
    for i in range(3):
        # We use the slider to scale the 'strength' of that specific mode
        mode_contribution = U[:, i] * S[i] * slider_vals[i]
        
        # Add to the base shape
        # Reshape mode_contribution if your barycenter is (N, 3) and U is flattened
        new_shape += mode_contribution.reshape(barycenter.shape)
        
    return new_shape


sliders = [1.2, -0.5, 0.1] # Tweak these to change the shape
new_deformed_shape = reconstruct_shape(
    U_soft_RL_centered, 
    S_soft_RL_centered, 
    barycenter_RL_stacked, 
    sliders
)


plot_latent_space_3d(S_soft_RL_centered, V_soft_RL_centered)



import pandas as pd
import numpy as np

# 1. Calculate the coordinates for each sample
# Scores = S * Vh (taking the first 3 rows of Vh and first 3 singular values)
# We transpose so each row is a sample and each column is a PC
latent_coords = (V_soft_RL_centered[:3, :].T * S_soft_RL_centered[:3])

# 2. Create a DataFrame
df_latent = pd.DataFrame(latent_coords, columns=['PC1', 'PC2', 'PC3'])

# 3. Add an index column for sample identification (optional)
df_latent['Sample_ID'] = range(len(df_latent))

# 4. Save to CSV
df_latent.to_csv('latent_space_data.csv', index=False)
print("Latent space coordinates saved to 'latent_space_data.csv'")
