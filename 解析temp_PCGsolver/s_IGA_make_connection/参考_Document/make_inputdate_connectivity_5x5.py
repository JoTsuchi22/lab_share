import numpy as np
import matplotlib.pyplot as plt
import math
import make_inputdate_function as mif

# difine filename
filename = "connectivity_5x5"

# patch info (ξ方向のコントロールポイント個数, η方向のコントロールポイント個数)
patch_info = np.array([5, 5])

# patch info 1行2列の時の例外処理
if np.array(patch_info.shape).shape[0] == 1:
    patch_info = np.array([[patch_info[0], patch_info[1]]])

xi_max = 0
eta_max = 0
for i in range(patch_info.shape[0]):
    if xi_max < patch_info[i][0]:
        xi_max = patch_info[i][0]
    if eta_max < patch_info[i][1]:
        eta_max = patch_info[i][1]

# パッチのコントロールポイント 座標，重み
patch = np.zeros((patch_info.shape[0], xi_max * eta_max, 3))
patch_bool = np.zeros((patch_info.shape[0], xi_max * eta_max, 3), dtype=bool)
for i in range(patch_info.shape[0]):
    patch_bool[i,:patch_info[i][0]*patch_info[i][1],:] = True

# 各パッチ入力
num = 0
patch[num,:patch_info[num][0]*patch_info[num][1],:] = np.array(
   [[1.0000000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.1666666666666667e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.5000000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.8333333333333333e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [2.0000000000000000e+00,	0.0000000000000000e+00,	1.0000000000000000e+00],
    [1.0000000000000000e+00,	2.6120387496374142e-01,	9.0236892706218264e-01],
    [1.1666666666666665e+00,	3.0473785412436494e-01,	9.0236892706218264e-01],
    [1.4999999999999998e+00,	3.9180581244561213e-01,	9.0236892706218264e-01],
    [1.8333333333333330e+00,	4.7887377076685916e-01,	9.0236892706218264e-01],
    [2.0000000000000000e+00,	5.2240774992748285e-01,	9.0236892706218264e-01],
    [7.3459080339013583e-01,	7.3459080339013583e-01,	8.3728154510363761e-01],
    [8.5702260395515850e-01,	8.5702260395515839e-01,	8.3728154510363761e-01],
    [1.1018862050852039e+00,	1.1018862050852039e+00,	8.3728154510363761e-01],
    [1.3467498062152490e+00,	1.3467498062152488e+00,	8.3728154510363761e-01],
    [1.4691816067802717e+00,	1.4691816067802717e+00,	8.3728154510363761e-01],
    [2.6120387496374148e-01,	1.0000000000000000e+00,	9.0236892706218252e-01],
    [3.0473785412436505e-01,	1.1666666666666667e+00,	9.0236892706218252e-01],
    [3.9180581244561224e-01,	1.5000000000000002e+00,	9.0236892706218252e-01],
    [4.7887377076685939e-01,	1.8333333333333333e+00,	9.0236892706218252e-01],
    [5.2240774992748296e-01,	2.0000000000000000e+00,	9.0236892706218252e-01],
    [0.0000000000000000e+00,	1.0000000000000000e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.1666666666666667e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.5000000000000000e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	1.8333333333333333e+00,	1.0000000000000000e+00],
    [0.0000000000000000e+00,	2.0000000000000000e+00,	1.0000000000000000e+00]])

# 重複ありの総コントロールポイント数
allpointnumber_overlapping = 0
for i in range(patch_info.shape[0]):
    allpointnumber_overlapping += patch_info[i][0] * patch_info[i][1]

globalpoint = np.zeros((allpointnumber_overlapping))
globalpoint_x = np.zeros((allpointnumber_overlapping))
globalpoint_y = np.zeros((allpointnumber_overlapping))
globalpoint_w = np.zeros((allpointnumber_overlapping))
globalpoint_bool = np.zeros((allpointnumber_overlapping), dtype = bool)

localpoint = np.zeros((patch_info.shape[0], xi_max * eta_max))
localpoint_bool = np.zeros((patch_info.shape[0], xi_max * eta_max), dtype=bool)
for i in range(patch_info.shape[0]):
    localpoint_bool[i,:patch_info[i][0]*patch_info[i][1]] = True

xi_and_eta_max = xi_max
if xi_max < eta_max:
    xi_and_eta_max = eta_max
A = np.zeros((patch_info.shape[0], 2, 2, xi_and_eta_max))
A_bool = np.zeros((patch_info.shape[0], 2, 2, xi_and_eta_max), dtype=bool)

# connect patch
xi = 0
eta = 1
positive = 0
negative = 1

patch_number = 0
globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool = \
    mif.connect_patch_arg_0boundary(patch_number, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool)

# 座標，パッチコネクティビティー，境界の辺をtxt出力
# [patch_number, xi_or_eta, 0_or_1(int)] (auto marge)
boundary_array_0 = np.array([[0, eta, 0]])

boundary_array_1 = np.array([[0, eta, 1]])

boundary_number = 0
mif.write_date_header(filename)
mif.write_date_localpoint(filename, patch_info, localpoint, localpoint_bool)
mif.write_date_globalpoint(filename, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w)
boundary_number = mif.write_boundary(filename, A, A_bool, boundary_array_0, boundary_number)
boundary_number = mif.write_boundary(filename, A, A_bool, boundary_array_1, boundary_number)


# 描写
color = np.array(["r", "g", "b", "c", "m", "y", "k"])
fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)

A = np.zeros((globalpoint[globalpoint_bool].shape[0],3))
for i in range(globalpoint[globalpoint_bool].shape[0]):
    A[i][0] = globalpoint[globalpoint_bool][i]
    A[i][1] = globalpoint_x[globalpoint_bool][i]
    A[i][2] = globalpoint_y[globalpoint_bool][i]
ax1.plot(A[:,1], A[:,2], c=color[0], marker="", linewidth=0.7)
for i in range(globalpoint[globalpoint_bool].shape[0]):
    ax1.text(A[i,1], A[i,2], str(int(A[i,0])), c=color[3], fontsize=6)


ax1.set_aspect('equal', adjustable='box')
ax1.set_axisbelow(True)
fig.set_figheight(9)
fig.set_figwidth(12)
ax1.grid()
ax1.set_xlim(-1, 3)
ax1.set_ylim(-1, 3)
plt.show()