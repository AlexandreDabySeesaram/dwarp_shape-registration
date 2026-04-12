import marimo

__generated_with = "0.22.0"
app = marimo.App()


@app.cell
def _():
    import meshio

    lungs = ["RL", "LL"]
    for lung in lungs:
        mesh = meshio.read(f"Meshes/Fine_morphed_sphere_{lung}.xml")
    
        meshio.write(
            f"Meshes/Fine_morphed_sphere_{lung}_Gmsh.msh", 
            mesh, 
            file_format="gmsh22", # Use "gmsh" for 4.1, "gmsh22" for 2.2
            binary=False
        )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
