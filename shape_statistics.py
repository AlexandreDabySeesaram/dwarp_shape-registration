import dolfin
import dolfin_warp as dwarp

mesh = dolfin.Mesh("Meshes/Coarse_sphere_RL.xml")

dwarp.warp(
        working_folder                              = "Results/barycenter", 
        working_basename                            = "Barycenter_mapping",
        images_folder                               = "None",
        images_basename                             = "None",
        mappings_folder                             = "Mappings",
        mappings_basename                           = "Mapping_coarse_sphere_RL",
        mesh                                        = mesh,
        warping_type                                = "barycenter"                                            , # registration
        kinematics_type                             = "full",
        nonlinearsolver                             = "gradient_descent",
        images_quadrature                           = 6,
        n_iter_max                                  = 100,
        # relax_type                                  = "gss",
        min_gradient_step                           = 1e-6,
        gradient_step                               = 1,
        write_VTU_files                             = True,
        write_VTU_files_with_preserved_connectivity = True,
        print_iterations                            = True,
        tol_dU                                      = 1e-3, 
        relax_n_iter_max                            = 30, 
        normalize_energies                          = False,
        gradient_type                               = "L2",
        continue_after_fail                         = True,
        inner_product_H1_weight                     = 1e-3)