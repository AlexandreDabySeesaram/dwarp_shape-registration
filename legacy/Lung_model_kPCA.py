import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk

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

#%% Fixed point

# On the full mapping (get also average homthety and average rigid body motion)

Mean_sphere_RL_array = mapping_sphere_RL_array.sum(axis=1)
Mean_sphere_LL_array = mapping_sphere_LL_array.sum(axis=1)

Mean_sphere_RL_array_centered = mapping_sphere_RL_array - Mean_sphere_RL_array[:,None]
Mean_sphere_LL_array_centered = mapping_sphere_LL_array - Mean_sphere_LL_array[:,None]

U_soft_sphere_RL_centered,S_soft_sphere_RL_centered,V_soft_sphere_RL_centered = np.linalg.svd(Mean_sphere_RL_array_centered)
U_soft_sphere_LL_centered,S_soft_sphere_LL_centered,V_soft_sphere_LL_centered = np.linalg.svd(Mean_sphere_LL_array_centered)

plt.figure()

plt.semilogy(S_soft_sphere_RL_centered[:-4], label = "Right lung")
plt.semilogy(S_soft_sphere_LL_centered[:-4], label = "Left lung")

plt.xlabel('Indexes of modes')
plt.ylabel('Singular values')
plt.legend()
plt.title("Spherical mappings w.r.t. mean shape")
plt.show()
plt.close()


x = np.linspace(0,S_soft_sphere_RL[:-4].shape[0],S_soft_sphere_RL[:-4].shape[0])
y_RL = S_soft_sphere_RL[:-4]
y_LL = S_soft_sphere_LL[:-4]
y_RL_full = S_sphere_RL[:-4]
y_LL_full = S_sphere_LL[:-4]
y_RL_centered = S_soft_sphere_RL_centered[:-4]
y_LL_centered = S_soft_sphere_LL_centered[:-4]

df = pd.DataFrame(np.stack((x,y_RL,y_LL,y_RL_full,y_LL_full,y_RL_centered,y_LL_centered)).T, columns=['N','SV_RL','SV_LL','SV_full_RL','SV_full_LL','SV_centered_RL','SV_centered_LL'])
df.to_csv('Results/Figure_Singular_values_mean.csv', index=False)


#%% Error SVD

N_rank = 0

error_vect = []

for N_rank in range(1,39):
    error_vect.append(100*np.linalg.norm(((U_soft_sphere_LL_centered[:,:N_rank]*S_soft_sphere_LL_centered[:N_rank]@V_soft_sphere_LL_centered[:N_rank,:]
) - Mean_sphere_LL_array_centered), axis=0)/np.linalg.norm(Mean_sphere_LL_array_centered, axis=0))



# plt.semilogy(np.stack(error_vect).T)
plt.semilogy(np.max(np.stack(error_vect), axis=1))

df = pd.DataFrame(np.stack((np.max(np.stack(error_vect), axis=1),np.linspace(1,38,38))).T, columns=['N','error_rel'])
df.to_csv('Results/SVD_erro.csv', index=False)

# Proj of first singular value in space

coef_project = U_soft_sphere_LL_centered.T@Mean_sphere_LL_array_centered

#%% kPCA 
def gaussian_kernel(X, sigma=1.0):
    """Compute the Gaussian (RBF) kernel matrix for a dataset X."""
    pairwise_sq_dists = np.sum(X**2, axis=1).reshape(-1, 1) + np.sum(X**2, axis=1) - 2 * np.dot(X, X.T)
    K = np.exp(-pairwise_sq_dists / (2 * sigma**2))
    return K

def center_kernel_matrix(K):
    """Center the kernel matrix K."""
    n = K.shape[0]
    one_n = np.ones((n, n)) / n
    K_centered = K - one_n @ K - K @ one_n + one_n @ K @ one_n
    return K_centered

def kpca(X, sigma=1.0, n_components=2):
    """Compute the top n_components of Kernel PCA using Gaussian kernel and SVD."""
    # Step 1: Compute the Gaussian kernel matrix
    K = gaussian_kernel(X, sigma=sigma)
    
    # Step 2: Center the Kernel matrix
    K_centered = center_kernel_matrix(K)
    
    # Step 3: Perform SVD on the centered kernel matrix
    U, S, Vt = np.linalg.svd(K_centered)
    
    # Step 4: Select the top n_components
    U_k = U[:, :n_components]  # top eigenvectors
    S_k = np.diag(S[:n_components])  # top singular values


    # Project the centered kernel matrix onto the top eigenvectors
    X_kpca = U_k @ S_k  # shape: (n_samples, n_components)
    
    return X_kpca


n_components = 3
sigma =1500  # Gamma for Gaussian kernel


Y_LL_centered_kpca = kpca(np.transpose(Mean_sphere_LL_array_centered), sigma=sigma, n_components=n_components)
print("Projected data shape:", Y_LL_centered_kpca.shape)
# %%











import plotly.graph_objects as go

# Sample data (replace this with your actual data)
# X = np.random.rand(3, 100)  # Example data
# colors[:, :3]

colors = plt.cm.viridis(np.linspace(0, 1, len(Y_LL_centered_kpca[:, 0])))                                        # Generate a color gradient


# Extract each axis from X
x_points = Y_LL_centered_kpca[:, 0]
y_points = Y_LL_centered_kpca[:, 1]
z_points = Y_LL_centered_kpca[:, 2]

# Create the 3D scatter plot
fig = go.Figure(data=[go.Scatter3d(
    x=x_points,
    y=y_points,
    z=z_points,
    mode='markers',
    marker=dict(
        size=4,  # Size of the markers
        color=colors[:,2],  # Color of the markers
        colorscale='viridis',  # Color scale
        showscale=True,  # Show color scale
    )
)])

# Set axes labels and title
fig.update_layout(
    scene=dict(
        xaxis_title='x',
        yaxis_title='y',
        zaxis_title='z',
    ),
    title='Physical Space'
)

# Show the plot
fig.show()
# %%


df = pd.DataFrame(np.stack((Y_LL_centered_kpca[:,0],Y_LL_centered_kpca[:,1],Y_LL_centered_kpca[:,2])).T, columns=['X1','X2','X3'])
df.to_csv('Results/kPCA_lungs.csv', index=False)



# %%
