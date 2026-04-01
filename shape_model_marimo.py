import marimo

__generated_with = "0.22.0"
app = marimo.App(width="columns")


@app.cell(column=0, hide_code=True)
def _(mo):
    mo.md(r"""
    This notebook proposes an interactive way to explore the "high" dimensional (more than 3 components) latent space composing the lungs shapes database

    ## Computing the latent space

    We first compute the latent space using an SVD on the centred mappings (mappings from which both the reduced kinematics and the shape barycenter have been deduced)
    """)
    return


@app.cell(hide_code=True)
def _(plotter):
    import vtk
    import numpy as np
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    import pyvista as pv
    import matplotlib.pyplot as plt
    job = 'configurations/barycenter_arm_config'



    def update_plot(mesh, mapping, initial_nodes):
        mesh.points = initial_nodes.copy()  
        mesh.points += mapping
        plotter.render()

    def vtu2numpy(vtu_name, return_output=False, Verbose = False):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(vtu_name + '.vtu')
        if Verbose:
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
        ax.set_title('3D Latent Space (Top 3 SVD modes)', fontsize=15)
        ax.set_xlabel(f'PC1 (Var: {S[0] ** 2 / np.sum(S ** 2):.1%})')
        ax.set_ylabel(f'PC2 (Var: {S[1] ** 2 / np.sum(S ** 2):.1%})')
        ax.set_zlabel(f'PC3 (Var: {S[2] ** 2 / np.sum(S ** 2):.1%})')
        fig.colorbar(scatter, ax=ax, label='PC1 Magnitude', shrink=0.5)
        plt.tight_layout()
        plt.savefig('latent_space_3d.png')
        plt.show()  

    N_patients = 40  
    Patients_Ids = list(range(0, N_patients + 1))  
    Lungs = ['RL']  
    mapping_RL = []
    indexes_mapping_RL = []
    mapping_LL = []
    indexes_mapping_LL = []
    reduced_mapping_RL = []
    reduced_mapping_LL = []
    reduced_mapping_LL = []  

    mapping_base_name = 'Mappings/Mapping_arm_Fine_mesh_morphed_sphere'
    for patient in Patients_Ids:
        for lung in Lungs:
            if lung == 'RL':  
                try:
                    mapping_RL.append(vtu2numpy(mapping_base_name + '_' + lung + '_' + str(patient).zfill(2)))
                    reduced_mapping_RL.append(vtu2numpy(mapping_base_name + '_reduced' + '_' + lung + '_' + str(patient).zfill(2)))  # Labeling
                    indexes_mapping_RL.append(patient)
                    print(f'Mapping {patient} accounted for')
                except:
                    print(f'Mapping {patient} not in database')
            else:
                try:
                    mapping_LL.append(vtu2numpy(mapping_base_name + '_' + lung) + str(patient).zfill(2))
                    reduced_mapping_LL.append(vtu2numpy(mapping_base_name + '_reduced' + '_' + lung) + str(patient).zfill(2))  # Add a colorbar to show the gradient of the first component
                    indexes_mapping_LL.append(patient)
                except:
                    print('Mapping not in database')
    barycenter_RL = vtu2numpy('Results/barycenter/barycenter_arm_Fine_mesh_morphed_sphere_ogdenciarletgeymonatneohookean_RL-frame=None_5130')
    barycenter_RL_stacked = barycenter_RL.flatten()[:, None]
    reduced_mapping_RL_stacked = [mapping.flatten()[:, None] for mapping in reduced_mapping_RL]
    mapping_RL_stacked = [mapping.flatten()[:, None] for mapping in mapping_RL]
    mapping_RL_array = np.hstack(mapping_RL_stacked)
    reduced_mapping_RL_array = np.hstack(reduced_mapping_RL_stacked)

    U_soft_RL_centered, S_soft_RL_centered, V_soft_RL_centered = np.linalg.svd(mapping_RL_array - reduced_mapping_RL_array - barycenter_RL_stacked)

    # plt.figure()
    # plt.semilogy(S_soft_RL_centered[:-4], label='Right lung')
    # plt.xlabel('Indexes of modes')
    # plt.ylabel('Singular values')
    # plt.legend()
    # plt.title('Spherical mappings full mapping')
    # plt.show()
    # plt.close()
    return (
        S_soft_RL_centered,
        V_soft_RL_centered,
        indexes_mapping_RL,
        lung,
        mapping_RL,
        mapping_base_name,
        pv,
        vtu2numpy,
    )


@app.cell(hide_code=True)
def _(S_soft_RL_centered):
    import plotly.graph_objects as go

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        y=S_soft_RL_centered[:-4],
        mode='lines+markers',
        name='Right lung',
        hovertemplate='Mode Index: %{x}<br>Singular Value: %{y:.2e}'
    ))

    fig.update_layout(
        title='Signular values decay - Morphed shpere mappings',
        xaxis_title='Indexes of modes',
        yaxis_title='Singular values',
        yaxis_type="log",  # This handles the 'semilogy' equivalent
        template='plotly_dark', # Optional: fits the marimo dark theme well
        legend=dict(yanchor="top", y=0.99, xanchor="right", x=0.99)
    )

    fig.show()
    return


@app.cell
def _(S_soft_RL_centered, V_soft_RL_centered, indexes_mapping_RL):

    import pandas as pd
    latent_coords = V_soft_RL_centered[:, :].T * S_soft_RL_centered[:]
    column_names = [f'PC{i+1}' for i in range(latent_coords.shape[1])]

    df_latent = pd.DataFrame(latent_coords, columns=column_names)
    df_latent.insert(0, 'Patient_ID', indexes_mapping_RL)

    df_latent.to_csv('high_dim_latent_space_data.csv', index=False)

    print("High dim latent space coordinates saved to 'latent_space_data.csv'")  
    return (df_latent,)


@app.cell
def _():
    return


@app.cell(column=1, hide_code=True)
def _(mo):
    mo.md(r"""
    ## Plot the latent space in a `ParallelCoordinates` plot

    Here we can easily to browse the "unfolded" latent space
    """)
    return


@app.cell
def _(df_pl):
    df_pl
    return


@app.cell(hide_code=True)
def _():
    import marimo as mo

    n_modes_slider = mo.ui.slider(start=3, stop=15, step=1, value=8, label="Number of modes")
    n_modes_slider
    return mo, n_modes_slider


@app.cell(hide_code=True)
def _(mo, n_modes_slider):
    available_pcs = [f"PC{i+1}" for i in range(n_modes_slider.value-1)]

    color_selector = mo.ui.dropdown(
        options=available_pcs, 
        value="PC1", 
        label="Colour by component"
    )

    color_selector
    return (color_selector,)


@app.cell(hide_code=True)
def _(color_selector, df_latent, mo, n_modes_slider):
    from wigglystuff import ParallelCoordinates
    import polars as pl

    # Convert pandas DataFrame to Polars 
    df_pl = pl.from_pandas(df_latent.iloc[:, :n_modes_slider.value]) 
    df_pl = pl.from_pandas(df_latent.iloc[:, :n_modes_slider.value]).rename({"Patient_ID": "uid"})
    # Interactive parallel coordinate plot
    parallel_plot = mo.ui.anywidget(
        ParallelCoordinates(
            df_pl, 
            height=500, 
            color_by=color_selector.value
        )
    )

    parallel_plot
    return df_pl, parallel_plot


@app.cell(hide_code=True)
def _(mo):
    zoom_slider = mo.ui.slider(start=0, stop=5, step=0.1, value=1, label="zoom level")
    zoom_slider
    return (zoom_slider,)


@app.cell(hide_code=True)
def _(mo):
    number_plots_slider = mo.ui.slider(start=1, stop=40, step=1, value=5, label="Number of pyvista plots")
    number_plots_slider
    return (number_plots_slider,)


@app.cell(hide_code=True)
def _(
    df_pl,
    lung,
    mapping_RL,
    mapping_base_name,
    mo,
    number_plots_slider,
    pv,
    vtu2numpy,
    zoom_slider,
):
    import tempfile
    def get_mesh_grid(indices):


        selected_indices = indices[:number_plots_slider.value]
        num_items = len(selected_indices)

        cols = 5
        rows = (num_items + cols - 1) // cols  

        plotter = pv.Plotter(
            shape=(rows, cols), 
            window_size=[900, 400 * rows], 
            off_screen=True, 
            notebook=True
        )
        for i, idx in enumerate(selected_indices):
            try:
                patient = df_pl["Patient_ID"][idx]
            except: 
                patient = df_pl["uid"][idx]
            disp, base_mesh = vtu2numpy(mapping_base_name + '_' + lung + '_' + str(patient),return_output=True)


            # Reshape to (N, 3) for PyVista
            disp = mapping_RL[idx].reshape(-1, 3)

            pv_mesh = pv.wrap(base_mesh) 

            current_mesh = pv_mesh.copy()
            current_mesh.points += disp

            # revert z-axis direction
            current_mesh.points[:, 2] *= -1

            # Add to subplot
            plotter.camera.zoom(zoom_slider.value)        
            plotter.subplot(i // cols, i % cols)
            plotter.add_mesh(current_mesh, color="lightblue", show_edges=False, smooth_shading=True)
            plotter.add_text(f"Sample Index: {idx}", font_size=10)
            plotter.view_isometric()

        # Return the interactive plotter window
        # return plotter.show(jupyter_backend='static', return_viewer=True)
        # return plotter.show(jupyter_backend='trame', return_viewer=True)

        # Output html from pyvista to display in marimo
        with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as tmp:
                plotter.export_html(tmp.name)

                with open(tmp.name, 'r') as f:
                    html_content = f.read()

        return mo.Html(f'<iframe srcdoc="{html_content.replace('"', '&quot;')}" width="100%" height="{400 * rows}px" style="border:none;"></iframe>')






    return (get_mesh_grid,)


@app.cell
def _():
    return


@app.cell(column=2)
def _(get_mesh_grid, parallel_plot):
    # import nest_asyncio
    # nest_asyncio.apply()


    # from pyvista.trame.jupyter import launch_server
    # await launch_server().ready

    get_mesh_grid(parallel_plot.widget.filtered_indices)
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
