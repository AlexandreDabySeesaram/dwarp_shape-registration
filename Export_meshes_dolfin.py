import meshio

# Load the .vtu file
mesh = meshio.read("Results/Mapping_sphere_PA37_RL-frame=0_007.vtu")

# Write to Dolfin XML
meshio.write("mesh.xml", mesh, file_format="dolfin-xml")

from dolfin import Mesh

mesh_2 = Mesh("mesh.xml")
