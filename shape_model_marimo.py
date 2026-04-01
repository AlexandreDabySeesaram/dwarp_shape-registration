import marimo

__generated_with = "0.22.0"
app = marimo.App(width="columns")


@app.cell(column=0)
def _(plotter):
    import vtk
    import numpy as np
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    import pyvista as pv
    import matplotlib.pyplot as plt
    job = 'configurations/barycenter_arm_config'
    # try:
    #     import tomllib
    # except ImportError:
    # with open(f'{job}.toml', 'rb') as f:
    #     config = tomllib.load(f)
      # Python 3.11+
    def update_plot(mesh, mapping, initial_nodes):
        mesh.points = initial_nodes.copy()  # Older versions
        mesh.points += mapping
        plotter.render()

    def vtu2numpy(vtu_name, return_output=False):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu_name + '.vtu')
        print(f"Name : {vtu_name + '.vtu'}")
        reader.Update()
        output = reader.GetOutput()
        displacement_vtk_array = output.GetPointData().GetArray('displacement')
        displacement_array = vtk_to_numpy(displacement_vtk_array)
        if return_output:
            return (displacement_array, output)
        else:
            return displacement_array

    def plot_latent_space_3d(S, Vh):
        """
        S: Singular values (S_soft_RL_centered)
        Vh: Right singular vectors (V_soft_RL_centered)
        """
        pc1 = S[0] * Vh[0, :]
        pc2 = S[1] * Vh[1, :]
        pc3 = S[2] * Vh[2, :]
        fig = plt.figure(figsize=(10, 7))
        ax = fig.add_subplot(111, projection='3d')
        scatter = ax.scatter(pc1, pc2, pc3, c=pc1, cmap='viridis', s=50, alpha=0.8)
        ax.set_title('3D Latent Space (Top 3 Principal Components)', fontsize=15)
        ax.set_xlabel(f'PC1 (Var: {S[0] ** 2 / np.sum(S ** 2):.1%})')
        ax.set_ylabel(f'PC2 (Var: {S[1] ** 2 / np.sum(S ** 2):.1%})')
        ax.set_zlabel(f'PC3 (Var: {S[2] ** 2 / np.sum(S ** 2):.1%})')
        fig.colorbar(scatter, ax=ax, label='PC1 Magnitude', shrink=0.5)
        plt.tight_layout()
        plt.savefig('latent_space_3d.png')
        plt.show()  # 1. Calculate the Latent Coordinates (Scores)
    N_patients = 40  # We take the first 3 components. 
    Patients_Ids = list(range(2, N_patients + 1))  # Coordinates = S_i * Vh_i
    Lungs = ['RL']  # This gives us the position of each 'sample' in the 3D reduced space.
    mapping_RL = []
    mapping_LL = []
    reduced_mapping_RL = []
    reduced_mapping_LL = []
    reduced_mapping_LL = []  # 2. Create the 3D Plot
    mapping_base_name = 'Mappings/Mapping_arm_Fine_mesh_morphed_sphere'
    for patient in Patients_Ids:
        for lung in Lungs:
            if lung == 'RL':  # Scatter plot
                try:
                    mapping_RL.append(vtu2numpy(mapping_base_name + '_' + lung + '_' + str(patient)))
                    reduced_mapping_RL.append(vtu2numpy(mapping_base_name + '_reduced' + '_' + lung + '_' + str(patient)))  # Labeling
                except:
                    print('Mapping not in database')
            else:
                try:
                    mapping_LL.append(vtu2numpy(mapping_base_name + '_' + lung) + str(patient))
                    reduced_mapping_LL.append(vtu2numpy(mapping_base_name + '_reduced' + '_' + lung) + str(patient))  # Add a colorbar to show the gradient of the first component
                except:
                    print('Mapping not in database')
    barycenter_RL = vtu2numpy('Results/barycenter/barycenter_arm_Fine_mesh_morphed_sphere_ogdenciarletgeymonatneohookean_RL-frame=None_5130')
    barycenter_RL_stacked = barycenter_RL.flatten()[:, None]
    reduced_mapping_RL_stacked = [mapping.flatten()[:, None] for mapping in reduced_mapping_RL]
    mapping_RL_stacked = [mapping.flatten()[:, None] for mapping in mapping_RL]
    mapping_RL_array = np.hstack(mapping_RL_stacked)
    reduced_mapping_RL_array = np.hstack(reduced_mapping_RL_stacked)

    U_soft_RL_centered, S_soft_RL_centered, V_soft_RL_centered = np.linalg.svd(mapping_RL_array - reduced_mapping_RL_array - barycenter_RL_stacked)
    # user input
    plt.figure()
    plt.semilogy(S_soft_RL_centered[:-4], label='Right lung')
    plt.xlabel('Indexes of modes')
    plt.ylabel('Singular values')
    plt.legend()
    plt.title('Spherical mappings full mapping')
    plt.show()
    plt.close()


    # plot_latent_space_3d(S_soft_RL_centered, V_soft_RL_centered)
    return S_soft_RL_centered, V_soft_RL_centered, mapping_RL_array


@app.cell
def _(S_soft_RL_centered, V_soft_RL_centered):

    import pandas as pd
    latent_coords = V_soft_RL_centered[:3, :].T * S_soft_RL_centered[:3]
    column_names = [f'PC{i+1}' for i in range(latent_coords.shape[1])]

    df_latent = pd.DataFrame(latent_coords, columns=column_names)

    df_latent['Sample_ID'] = range(len(df_latent))
    df_latent.to_csv('high_dim_latent_space_data.csv', index=False)

    print("High dim latent space coordinates saved to 'latent_space_data.csv'")  
    return


@app.cell
def _(mapping_RL_array):
    mapping_RL_array.shape
    return


@app.cell
def _():
    return


@app.cell(column=1)
def _():
    return


if __name__ == "__main__":
    app.run()
