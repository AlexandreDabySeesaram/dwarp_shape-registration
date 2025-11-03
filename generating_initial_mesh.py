import meshio
import numpy as np


lungs       = ["RL", "LL"]
fineness    = "fine" 
coef_relax  = 0.5


match fineness:
    case "coarse":
        for lung in lungs:
            mesh = meshio.read("Results/Mapping_sphere_PA37_"+lung+"_000.vtu")

            morphed_mesh = meshio.Mesh(
                points=mesh.points + coef_relax*mesh.point_data.get("displacement"),
                cells=mesh.cells
            )

            # Write to Dolfin XML
            meshio.write(
                "Meshes/Coarse_morphed_sphere_"+lung+".xml", morphed_mesh, file_format="dolfin-xml")
            print(f"* Finished writing {fineness} {lung}.")
    case "fine":
        basename = "../rMartin_results/Lungs_mapping"
        for lung in lungs:
            mesh = meshio.read(basename+"/Mapping_fine_mesh_sphere_PA37_"+lung+"_000.vtu")

            morphed_mesh = meshio.Mesh(
                points=mesh.points + coef_relax*mesh.point_data.get("displacement"),
                cells=mesh.cells
            )

            # Write to Dolfin XML
            meshio.write(
                "Meshes/Fine_morphed_sphere_"+lung+".xml", morphed_mesh, file_format="dolfin-xml")
            print(f"* Finished writing {fineness} {lung}.")


