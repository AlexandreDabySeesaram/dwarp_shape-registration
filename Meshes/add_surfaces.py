import gmsh
import sys

gmsh.initialize()
gmsh.open("Fine_morphed_sphere_RL_Gmsh.msh") 
# 1. Get all 3D entities (volumes)
volumes = gmsh.model.getEntities(3)


# 1. Get all 3D volumes
volumes = gmsh.model.getEntities(3)
vol_tags = [v[1] for v in volumes]

# 2. Get the boundary of the volumes
# This returns the 2D surface entities that bound the 3D volumes
boundary = gmsh.model.getBoundary(volumes, combined=True, oriented=False, recursive=False)

# # # 3. Create a Physical Group for these surfaces
# # # This tells Gmsh to include these 2D elements in the export
# surf_tags = [b[1] for b in boundary]
# gmsh.model.addPhysicalGroup(2, surf_tags, 111)


# # Assign 111 to the 2D boundary elements
# gmsh.model.addPhysicalGroup(2, surf_tags, 111)

# # Assign an ID (e.g., 1) to the 3D bulk elements so they aren't 0 either
# gmsh.model.addPhysicalGroup(3, vol_tags, 1)


# 4. Crucial step: If the mesh doesn't have 2D elements yet,
# we force the creation of boundary elements without 'createGeometry'
gmsh.model.mesh.createTopology()
# We skip classifySurfaces and createGeometry to avoid the error in your screenshot


# 2. Get the boundary of the volumes
# This returns the 2D surface entities that bound the 3D volumes
boundary = gmsh.model.getBoundary(volumes, combined=True, oriented=False, recursive=False)

# # 3. Create a Physical Group for these surfaces
# # This tells Gmsh to include these 2D elements in the export
surf_tags = [b[1] for b in boundary]
gmsh.model.addPhysicalGroup(2, surf_tags, 111)


# Assign an ID (e.g., 1) to the 3D bulk elements so they aren't 0 either
gmsh.model.addPhysicalGroup(3, vol_tags, 1)


# 5. Save the file
# This will export the original 3D elements + the new 2D boundary elements
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("new_mesh_with_skin.msh")

gmsh.finalize()
