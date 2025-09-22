import meshio
import glob
import dolfin
import dolfin_mech as dmech
mesh_RL = dolfin.Mesh("Meshes/Fine_sphere_RL.xml")  

def interpolate_vtu_to_dolfin(vtu_path, V, target_mesh):
    # Read VTU file
    vtu_data = meshio.read(vtu_path)
    points = vtu_data.points
    displacements = vtu_data.point_data["displacement"]  # Adjust field name

    # Create a temporary Function for source data
    source_V = dolfin.VectorFunctionSpace(target_mesh, "Lagrange", 1)
    u_source = dolfin.Function(source_V)
    u_source.vector()[:] = displacements.flatten()

    # Interpolate to target space
    u_target = dolfin.Function(V,name="displacement")
    u_target.interpolate(u_source)
    return u_target



mappings_folder = "Mappings"
mappings_basename = "Mapping_fine_sphere_RL"
mappings_ext = ".vtu"
tag_interpolation = "interpolated"

initial_names_list = glob.glob(mappings_folder+"/"+mappings_basename+"_??"+mappings_ext)

V = dolfin.VectorFunctionSpace(mesh_RL, "Lagrange", 1)  # P1 for displacement

for mapping in initial_names_list:

    u_interpolated = interpolate_vtu_to_dolfin(mapping, V, mesh_RL)
    new_name = mapping.replace("sphere_", "sphere_"+tag_interpolation+"_")
    new_name = new_name.replace(".vtu", "")
    dmech.write_VTU_file(
        filebasename = new_name,
        function = u_interpolated,
        preserve_connectivity = True)