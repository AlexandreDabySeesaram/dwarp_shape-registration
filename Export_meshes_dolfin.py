import meshio

fineness = "fine"

match fineness:
    case "coarse":
        # Load the .vtu file
        mesh_RL = meshio.read("Results/Mapping_sphere_PA37_RL-frame=0_007.vtu")
        mesh_LL = meshio.read("Results/Mapping_sphere_PA37_LL-frame=0_007.vtu")

        # Write to Dolfin XML
        meshio.write("Meshes/Coarse_sphere_RL.xml", mesh_RL, file_format="dolfin-xml")
        meshio.write("Meshes/Coarse_sphere_LL.xml", mesh_LL, file_format="dolfin-xml")


        import os
        import shutil

        # make sure destination folder exists
        os.makedirs("Mappings", exist_ok=True)

        for pid in range(2, 41):  # patient IDs 2 → 40
            pid_str = f"{pid:02d}"  # zero-padded, e.g. 02, 03 ...
            src = f"Results/Mapping_sphere_PA{pid}_RL_000.vtu"
            dst = f"Mappings/Mapping_coarse_sphere_RL_{pid_str}.vtu"

            if os.path.exists(src):
                shutil.copy(src, dst)
                print(f"Copied {src} -> {dst}")
            else:
                print(f"⚠️ File not found: {src}")
    case "fine":
        basename = "../rMartin_results/Lungs_mapping"
        # Load the .vtu file
        mesh_RL = meshio.read(basename+"/Mapping_fine_mesh_sphere_PA37_RL-frame=0_007.vtu")
        mesh_LL = meshio.read(basename+"/Mapping_fine_mesh_sphere_PA37_LL-frame=0_007.vtu")

        # Write to Dolfin XML
        meshio.write("Meshes/Fine_sphere_RL.xml", mesh_RL, file_format="dolfin-xml")
        meshio.write("Meshes/Fine_sphere_LL.xml", mesh_LL, file_format="dolfin-xml")


        import os
        import shutil

        # make sure destination folder exists
        os.makedirs("Mappings", exist_ok=True)

        for pid in range(2, 41):  # patient IDs 2 → 40
            pid_str = f"{pid:02d}"  # zero-padded, e.g. 02, 03 ...
            src_RL = basename+f"/Mapping_fine_mesh_sphere_PA{pid}_RL_000.vtu"
            dst_RL = f"Mappings/Mapping_fine_sphere_RL_{pid_str}.vtu"
            src_LL = basename+f"/Mapping_fine_mesh_sphere_PA{pid}_LL_000.vtu"
            dst_LL = f"Mappings/Mapping_fine_sphere_LL_{pid_str}.vtu"



            if os.path.exists(src_RL):
                shutil.copy(src_RL, dst_RL)
                print(f"Copied {src_RL} -> {dst_RL}")
            else:
                print(f"⚠️ File not found: {src_RL}")
            if os.path.exists(src_LL):
                shutil.copy(src_LL, dst_LL)
                print(f"Copied {src_LL} -> {dst_LL}")
            else:
                print(f"⚠️ File not found: {src_LL}")