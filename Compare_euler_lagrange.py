import numpy as np
import matplotlib.pyplot as plt


# grad_euler_name = "/Users/daby/Documents/Code/Shape_registration/res_euler.dat"
# grad_lagrange_name = "res_lagrange.dat"

# grad_euler = np.loadtxt(grad_euler_name)
# grad_lagrange = np.loadtxt(grad_lagrange_name)

# Delta_grad = np.sqrt(np.sum((grad_euler-grad_lagrange)**2,axis=0))/np.sqrt(np.sum((grad_lagrange)**2,axis=0))

# plt.figure()

# plt.semilogy(Delta_grad)
# plt.ylabel("Relative sobolev gradient error (log)")
# plt.xlabel("Iterations")
# plt.show()




shape_deriv_euler_name = "/Users/daby/Documents/Code/Shape_registration/form_euler.npy"
shape_deriv_lagrange_name = "form_lagrange_noF.npy"



shape_deriv_euler = np.load(shape_deriv_euler_name)
shape_deriv_lagrange = np.load(shape_deriv_lagrange_name)



Delta_shape_deriv = np.linalg.norm(shape_deriv_euler-shape_deriv_lagrange, axis=0)/np.linalg.norm(shape_deriv_euler, axis=0)


plt.figure()
plt.semilogy(Delta_shape_deriv)
plt.ylabel("Relative gateau derivative error (log)")
plt.xlabel("Iterations")
plt.show()





inner_product_euler_name = "/Users/daby/Documents/Code/Shape_registration/form_bi_euler.npy"
inner_product_euler_name_1 = "/Users/daby/Documents/Code/Shape_registration/form_bi_euler_1.npy"
inner_product_euler_name_2 = "/Users/daby/Documents/Code/Shape_registration/form_bi_euler_2.npy"
# inner_product_lagrange_name = "form_bi_lagrange.npy"
inner_product_lagrange_name_noF = "form_bi_lagrange_noF.npy"
inner_product_lagrange_name_noF_1 = "form_bi_lagrange_1_noF.npy"
inner_product_lagrange_name_noF_2 = "form_bi_lagrange_2_noF.npy"



inner_product_euler = np.load(inner_product_euler_name)
inner_product_euler_1 = np.load(inner_product_euler_name_1)
inner_product_euler_2 = np.load(inner_product_euler_name_2)
# inner_product_lagrange = np.load(inner_product_lagrange_name)
inner_product_lagrange_noF = np.load(inner_product_lagrange_name_noF)
inner_product_lagrange_noF_1 = np.load(inner_product_lagrange_name_noF_1)
inner_product_lagrange_noF_2 = np.load(inner_product_lagrange_name_noF_2)

# Delta_inner_product = np.linalg.norm(inner_product_euler-inner_product_lagrange, axis=(0,1))/np.linalg.norm(inner_product_euler, axis=(0,1))
Delta_inner_product_noF = np.linalg.norm(inner_product_euler-inner_product_lagrange_noF, axis=(0,1))/np.linalg.norm(inner_product_euler, axis=(0,1))
Delta_inner_product_noF_1 = np.linalg.norm(inner_product_euler_1-inner_product_lagrange_noF_1, axis=(0,1))/np.linalg.norm(inner_product_euler_1, axis=(0,1))
Delta_inner_product_noF_2 = np.linalg.norm(inner_product_euler_2-inner_product_lagrange_noF_2, axis=(0,1))/np.linalg.norm(inner_product_euler_2, axis=(0,1))
plt.figure()

# plt.semilogy(Delta_inner_product)
# plt.ylabel("Relative inner product error (log)")
# plt.xlabel("Iterations")
# plt.show()

plt.semilogy(Delta_inner_product_noF)
plt.ylabel("Relative inner product error (log)")
plt.xlabel("Iterations")
plt.show()

plt.semilogy(Delta_inner_product_noF_1)
plt.ylabel("Relative inner product error semi H1 (log)")
plt.xlabel("Iterations")
plt.show()

plt.semilogy(Delta_inner_product_noF_2)
plt.ylabel("Relative inner product L2 error (log)")
plt.xlabel("Iterations")
plt.show()



inner_product_F_name = "F_inner_lagrange.npy"
inner_product_F = np.load(inner_product_F_name)


plt.plot(inner_product_F)
plt.ylabel("norm F_inv")
plt.xlabel("Iterations")
plt.show()