import dolfin_mech as dmech
import dolfin_warp as dwarp
import NL_eneregies










problem = dwarp.FullKinematicsWarpingProblem(
    mesh=mesh,
    mesh_folder=mesh_folder,
    mesh_basename=mesh_basename,
    U_degree=mesh_degree,
    silent=silent)






barycentric_energy = NL_eneregies.BarycentricEnergy(
    problem=problem,
    images_series=images_series,
    quadrature_degree=images_quadrature,
    w=image_w,
    w_char_func=images_char_func,
    im_is_cone=images_is_cone,
    static_scaling=images_static_scaling,
    dynamic_scaling=images_dynamic_scaling)
problem.add_image_energy(warped_image_energy)








solver = dwarp.GradientDescentSolver(
    problem=problem,
    parameters={
        "relax_type"                : relax_type                ,
        "working_folder"            : working_folder            ,
        "working_basename"          : working_basename          ,
        "write_iterations"          : print_iterations          , 
        "min_gradient_step"         : min_gradient_step         , 
        "gradient_step"             : gradient_step             ,
        "n_iter_max"                : n_iter_max                ,
        "relax_n_iter_max"          : relax_n_iter_max          ,
        "tol_dU"                    : tol_dU                    ,
        "gradient_type"             : gradient_type             ,
        "inner_product_H1_weight"   : inner_product_H1_weight   ,
        })