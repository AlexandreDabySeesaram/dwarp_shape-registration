# import meshio
# import numpy as np
# import dolfin 


# #%% Reading meshes

# mesh_io = meshio.read("Results/Mapping_sphere_PA2_RL_000.vtu")    
# meshio.write("xml_sphere.xml", mesh_io)
# dolfin_mesh = dolfin.Mesh("xml_sphere.xml")


# #%% Extracting the operators



# # Define function space
# V = dolfin.VectorFunctionSpace(dolfin_mesh, 'P', 1)  # Using linear Lagrange vector elements

# # Define trial and test functions
# u = dolfin.TrialFunction(V)
# v = dolfin.TestFunction(V)

# # Define the bilinear form for the stiffness matrix
# a = dolfin.inner(dolfin.sym(dolfin.grad(u)), dolfin.sym(dolfin.grad(v))) * dolfin.dx

# # Assemble the stiffness matrix
# K = dolfin.assemble(a)

# # Optionally, you can convert the matrix to a dense format for inspection
# K_dense = K.array()

# # Print the shape of the stiffness matrix
# print("Shape of the stiffness matrix:", K_dense.shape)

# # Optionally, save the stiffness matrix to a file
# np.save("stiffness_matrix.npy", K_dense)


# #%% Load mappings and compute energies

# soft_mapping_flat = np.load("soft_mappings_flat.npy")
# K_dense = np.load("stiffness_matrix.npy")

import torch


K = torch.tensor(K_dense, dtype=torch.float32)
mappings = torch.tensor(soft_mapping_flat, dtype=torch.float32)
N = mappings.shape[0]

# E = soft_mapping_flat[:,0].T@A_dense@soft_mapping_flat[:,0]


def strain_energy(K,U):
    return 0.5 * torch.stack([torch.dot(u, torch.matmul(K, u)) for u in U.T])

def loss_function(U_bar, mappings, K):

    return torch.sum(strain_energy(K,U_bar - mappings))

#%% Optimisation


# Dimensions

# Initialize the vector with random values
U_bar = torch.randn(N, 1, requires_grad=True)

# Set up the optimizer
learning_rate = 0.01
optimizer = torch.optim.Adam([U_bar], lr=learning_rate)

Loss_vect = []

# Optimization loop
num_iterations = 1000
for i in range(num_iterations):
    # Zero the gradients
    optimizer.zero_grad()

    # Compute the loss
    loss = loss_function(U_bar, mappings, K)

    # Compute gradients
    loss.backward()

    # Update the vector
    optimizer.step()

    Loss_vect.append(loss.item())

    # Print progress
    if (i + 1) % 10 == 0:
        print(f'Iteration {i + 1}, Loss: {loss.item()}')

# Print the optimized vector
print("Optimized Vector:")
print(U_bar)

np.save("Average_mapping.npy", U_bar)
