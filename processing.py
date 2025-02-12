import dolfin_warp as dwarp
import dolfin



def reduced_kinematics(
        result_folder               : str               = None,
        filebasename                : str               = None,
        image_folder                : str               = None,
        image_name                  : str               = None,
        mesh                        : dolfin.Mesh       = None,
        reduced_kinematics_model    : str               = None,
        tol_dU                      : float             = 1e-3,
        n_iter_max                  : int               = 100,
        images_quadrature           : int               = 6,
        ):

    dwarp.warp(
            warping_type                                = "registration"                                            , # registration
            working_folder                              = result_folder,
            images_char_func                            = False,
            working_basename                            = filebasename+"_reduced",
            images_folder                               = image_folder,
            images_basename                             = image_name,
            images_ext                                  = "vti",
            mesh                                        = mesh,
            kinematics_type                             = "reduced",
            reduced_kinematics_model                    = reduced_kinematics_model,
            images_quadrature                           = 6,
            n_iter_max                                  = 100,
            nonlinearsolver                             = "gradient_descent",
            min_gradient_step                           = 1e-6,
            gradient_step                               = 1,
            continue_after_fail                         = 1,
            write_VTU_files                             = True,
            write_VTU_files_with_preserved_connectivity = True,
            initialize_reduced_U_from_file              = False,
            print_iterations                            = True,
            tol_dU                                      = 1e-3, 
            register_ref_frame                          = True,
            relax_n_iter_max                            = 30, 
            normalize_energies                          = False,
            gradient_type                               = "Sobolev",
            inner_product_H1_weight                     = 1e-3) #L2


def full_kinematics(
        result_folder               : str               = None,
        filebasename                : str               = None,
        image_folder                : str               = None,
        image_name                  : str               = None,
        mesh                        : dolfin.Mesh       = None,
        tol_dU                      : float             = 1e-3,
        n_iter_max                  : int               = 100,
        images_quadrature           : int               = 6,
        ):

    dwarp.warp(
            warping_type                                = "registration"                                            , # registration
            working_folder                              = result_folder,
            images_char_func                            = False,
            working_basename                            = filebasename,
            images_folder                               = image_folder,
            images_basename                             = image_name,
            images_ext                                  = "vti",
            mesh                                        = mesh,
            kinematics_type                             = "full",
            nonlinearsolver                             = "gradient_descent",
            images_quadrature                           = 6,
            n_iter_max                                  = 100,
            # relax_type                                  = "gss",
            min_gradient_step                           = 1e-6,
            gradient_step                               = 1,
            continue_after_fail                         = 1,
            write_VTU_files                             = True,
            write_VTU_files_with_preserved_connectivity = True,
            initialize_reduced_U_from_file              = False,
            print_iterations                            = True,
            tol_dU                                      = 1e-3, 
            register_ref_frame                          = True,
            relax_n_iter_max                            = 30, 
            normalize_energies                          = False,
            gradient_type                               = "Sobolev",
            inner_product_H1_weight                     = 1e-3,
            initialize_U_from_file                      = True,
            initialize_U_folder                         = result_folder, 
            initialize_U_basename                       = filebasename+"_reduced",
            initialize_U_ext                            = "vtu",
            initialize_U_array_name                     = "displacement",
            initialize_U_method                         = "interpolation") # dofs_transfer, interpolation, projection) #L2