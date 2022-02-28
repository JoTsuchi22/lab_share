import numpy as np
import matplotlib.pyplot as plt
import math

def connect_patch_arg_0boundary(patch_number, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool):
    length = int(globalpoint[globalpoint_bool].shape[0])
    a_g = length
    a_l_g = 0
    a_l_l = 0
    patch_x = np.reshape(patch[patch_number,:,0][patch_bool[patch_number,:,0]], [patch_info[patch_number][1], patch_info[patch_number][0]]).T
    patch_y = np.reshape(patch[patch_number,:,1][patch_bool[patch_number,:,1]], [patch_info[patch_number][1], patch_info[patch_number][0]]).T
    patch_w = np.reshape(patch[patch_number,:,2][patch_bool[patch_number,:,2]], [patch_info[patch_number][1], patch_info[patch_number][0]]).T
    for j in range(patch_info[patch_number][1]):
        for i in range(patch_info[patch_number][0]):
            localpoint[patch_number,a_l_g] = length + a_l_l
            localpoint_bool[patch_number,a_l_g] = True
            globalpoint[a_g] = a_g
            globalpoint_bool[a_g] = True
            globalpoint_x[a_g] = patch_x[i][j]
            globalpoint_y[a_g] = patch_y[i][j]
            globalpoint_w[a_g] = patch_w[i][j]
            a_g = a_g + 1
            a_l_g = a_l_g + 1
            a_l_l = a_l_l + 1
    a_l_g = 0
    for j in range(patch_info[patch_number][1]):
        for i in range(patch_info[patch_number][0]):
            if i == 0:
                A[patch_number,0,0,j] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,0,0,j] = True
            if i == patch_info[patch_number][0] - 1:
                A[patch_number,0,1,j] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,0,1,j] = True
            if j == 0:
                A[patch_number,1,0,i] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,1,0,i] = True
            if j == patch_info[patch_number][1] - 1:
                A[patch_number,1,1,i] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,1,1,i] = True
            a_l_g = a_l_g + 1
    return globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool


def connect_patch_arg_1boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool):
    xi = 0
    eta = 1
    positive = 0
    negative = 1

    if connection_vec[1] == xi:
        end_0 = patch_info[patch_number][0] - 1
    if connection_vec[1] == eta:
        end_0 = patch_info[patch_number][1] - 1
    if connection_vec[2] == 0:
        boundary = 0
    if connection_vec[2] == 1:
        boundary = end_0
    if connection_vec[4] == xi:
        axis = xi
        end_1 = patch_info[connection_vec[3]][0] - 1
    if connection_vec[4] == eta:
        axis = eta
        end_1 = patch_info[connection_vec[3]][1] - 1
    
    length = int(globalpoint[globalpoint_bool].shape[0])
    a_g = length
    a_l_g = 0
    a_l_l = 0
    patch_x = np.reshape(patch[patch_number,:,0][patch_bool[patch_number,:,0]], [patch_info[patch_number][1], patch_info[patch_number][0]]).T
    patch_y = np.reshape(patch[patch_number,:,1][patch_bool[patch_number,:,1]], [patch_info[patch_number][1], patch_info[patch_number][0]]).T
    patch_w = np.reshape(patch[patch_number,:,2][patch_bool[patch_number,:,2]], [patch_info[patch_number][1], patch_info[patch_number][0]]).T
    for j in range(patch_info[patch_number][1]):
        for i in range(patch_info[patch_number][0]):
            if connection_vec[1] == xi:
                point = i
                point_var = j
            if connection_vec[1] == eta:
                point = j
                point_var = i
            if connection_vec[6] == positive:
                direction = point_var
            if connection_vec[6] == negative:
                direction = end_1 - point_var
            if point == boundary:
                localpoint[patch_number,a_l_g] = A[connection_vec[3],axis,connection_vec[5],direction]
                localpoint_bool[patch_number,a_l_g] = True
                a_l_g = a_l_g + 1
            if point != boundary:
                localpoint[patch_number,a_l_g] = length + a_l_l
                localpoint_bool[patch_number,a_l_g] = True
                globalpoint[a_g] = a_g
                globalpoint_bool[a_g] = True
                globalpoint_x[a_g] = patch_x[i][j]
                globalpoint_y[a_g] = patch_y[i][j]
                globalpoint_w[a_g] = patch_w[i][j]
                a_g = a_g + 1
                a_l_g = a_l_g + 1
                a_l_l = a_l_l + 1
    a_l_g = 0
    for j in range(patch_info[patch_number][1]):
        for i in range(patch_info[patch_number][0]):
            if i == 0:
                A[patch_number,0,0,j] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,0,0,j] = True
            if i == patch_info[patch_number][0] - 1:
                A[patch_number,0,1,j] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,0,1,j] = True
            if j == 0:
                A[patch_number,1,0,i] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,1,0,i] = True
            if j == patch_info[patch_number][1] - 1:
                A[patch_number,1,1,i] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,1,1,i] = True
            a_l_g = a_l_g + 1
    return globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool


def connect_patch_arg_2boundary(patch_number, connection_vec, patch_info, patch, patch_bool, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool):
    xi = 0
    eta = 1
    positive = 0
    negative = 1

    if connection_vec[0][1] == xi:
        end_0 = patch_info[patch_number][0] - 1
    if connection_vec[0][1] == eta:
        end_0 = patch_info[patch_number][1] - 1
    if connection_vec[0][2] == 0:
        boundary_0 = 0
    if connection_vec[0][2] == 1:
        boundary_0 = end_0
    if connection_vec[0][4] == xi:
        axis_0 = xi
        end_1 = patch_info[connection_vec[0][3]][0] - 1
    if connection_vec[0][4] == eta:
        axis_0 = eta
        end_1 = patch_info[connection_vec[0][3]][1] - 1
    
    if connection_vec[1][1] == xi:
        end_2 = patch_info[patch_number][0] - 1
    if connection_vec[1][1] == eta:
        end_2 = patch_info[patch_number][1] - 1
    if connection_vec[1][2] == 0:
        boundary_1 = 0
    if connection_vec[1][2] == 1:
        boundary_1 = end_2
    if connection_vec[1][4] == xi:
        axis_1 = xi
        end_3 = patch_info[connection_vec[1][3]][0] - 1
    if connection_vec[1][4] == eta:
        axis_1 = eta
        end_3 = patch_info[connection_vec[1][3]][1] - 1

    length = int(globalpoint[globalpoint_bool].shape[0])
    a_g = length
    a_l_g = 0
    a_l_l = 0
    patch_x = np.reshape(patch[patch_number,:,0][patch_bool[patch_number,:,0]], [patch_info[patch_number][1], patch_info[patch_number][0]]).T
    patch_y = np.reshape(patch[patch_number,:,1][patch_bool[patch_number,:,1]], [patch_info[patch_number][1], patch_info[patch_number][0]]).T
    patch_w = np.reshape(patch[patch_number,:,2][patch_bool[patch_number,:,2]], [patch_info[patch_number][1], patch_info[patch_number][0]]).T
    for j in range(patch_info[patch_number][1]):
        for i in range(patch_info[patch_number][0]):
            if connection_vec[0][1] == xi:
                point_0 = i
                point_var_0 = j
            if connection_vec[0][1] == eta:
                point_0 = j
                point_var_0 = i
            if connection_vec[0][6] == positive:
                direction_0 = point_var_0
            if connection_vec[0][6] == negative:
                direction_0 = end_1 - point_var_0
            if connection_vec[1][1] == xi:
                point_1 = i
                point_var_1 = j
            if connection_vec[1][1] == eta:
                point_1 = j
                point_var_1 = i
            if connection_vec[1][6] == positive:
                direction_1 = point_var_1
            if connection_vec[1][6] == negative:
                direction_1 = end_3 - point_var_1
            if point_0 == boundary_0 and point_1 == boundary_1:
                localpoint[patch_number,a_l_g] = A[connection_vec[0][3],axis_0,connection_vec[0][5],direction_0]
                localpoint_bool[patch_number,a_l_g] = True
                a_l_g = a_l_g + 1
            if point_0 == boundary_0 and point_1 != boundary_1:
                localpoint[patch_number,a_l_g] = A[connection_vec[0][3],axis_0,connection_vec[0][5],direction_0]
                localpoint_bool[patch_number,a_l_g] = True
                a_l_g = a_l_g + 1
            if point_0 != boundary_0 and point_1 == boundary_1:
                localpoint[patch_number,a_l_g] = A[connection_vec[1][3],axis_1,connection_vec[1][5],direction_1]
                localpoint_bool[patch_number,a_l_g] = True
                a_l_g = a_l_g + 1
            if point_0 != boundary_0 and point_1 != boundary_1:
                localpoint[patch_number,a_l_g] = length + a_l_l
                localpoint_bool[patch_number,a_l_g] = True
                globalpoint[a_g] = a_g
                globalpoint_bool[a_g] = True
                globalpoint_x[a_g] = patch_x[i][j]
                globalpoint_y[a_g] = patch_y[i][j]
                globalpoint_w[a_g] = patch_w[i][j]
                a_g = a_g + 1
                a_l_g = a_l_g + 1
                a_l_l = a_l_l + 1
    a_l_g = 0
    for j in range(patch_info[patch_number][1]):
        for i in range(patch_info[patch_number][0]):
            if i == 0:
                A[patch_number,0,0,j] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,0,0,j] = True
            if i == patch_info[patch_number][0] - 1:
                A[patch_number,0,1,j] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,0,1,j] = True
            if j == 0:
                A[patch_number,1,0,i] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,1,0,i] = True
            if j == patch_info[patch_number][1] - 1:
                A[patch_number,1,1,i] = localpoint[patch_number,a_l_g]
                A_bool[patch_number,1,1,i] = True
            a_l_g = a_l_g + 1
    return globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w, localpoint, localpoint_bool, A, A_bool


# 3boundary 必要?? 必要だったら自分で作ってください，多分5分くらいでできる
# def connect_patch_arg_3boundary(patch_number, patch_info, patch):

# テキストデータ書き込み
def write_date_header(filename):
    filename_txt = filename + ".txt"
    f = open(filename_txt, 'w')
    f.write('patch connectivity')
    f.write('\n')


def write_date_localpoint(filename, patch_info, localpoint, localpoint_bool):
    filename_txt = filename + ".txt"
    f = open(filename_txt, 'a')
    for i in range(patch_info.shape[0]):
        localpoint_var = localpoint[i,:][localpoint_bool[i,:]]
        for j in range(localpoint_var.shape[0]):
            if j == 0:
                f.write(str(int(localpoint_var[j])))
            else:
                f.write('   ')
                f.write(str(int(localpoint_var[j])))
        f.write('\n')


def write_date_globalpoint(filename, globalpoint, globalpoint_bool, globalpoint_x, globalpoint_y, globalpoint_w):
    filename_txt = filename + ".txt"
    f = open(filename_txt, 'a')
    f.write('\n')
    f.write('control point infomation')
    f.write('\n')
    A = np.zeros((globalpoint[globalpoint_bool].shape[0], 4))
    for i in range(globalpoint[globalpoint_bool].shape[0]):
        A[i][0] = globalpoint[globalpoint_bool][i]
        A[i][1] = globalpoint_x[globalpoint_bool][i]
        A[i][2] = globalpoint_y[globalpoint_bool][i]
        A[i][3] = globalpoint_w[globalpoint_bool][i]
        f.write(str(int(A[i][0])))
        f.write('   ')
        f.write('   ')
        f.write(str('{:.21e}'.format(A[i][1])))
        f.write('   ')
        f.write('   ')
        f.write(str('{:.21e}'.format(A[i][2])))
        f.write('   ')
        f.write('   ')
        f.write(str('{:.21e}'.format(A[i][3])))
        f.write('\n')
    f.write('\n')


def write_boundary(filename, A, A_bool, boundary_array, boundary_number):
    filename_txt = filename + ".txt"
    f = open(filename_txt, 'a')
    f.write('boundary ' + str(boundary_number))
    f.write('\n')
    if (np.array(boundary_array.shape)).shape[0] == 1:
        B = A[boundary_array[0],boundary_array[1],boundary_array[2],:][A_bool[boundary_array[0],boundary_array[1],boundary_array[2],:]]
        B = np.sort(B)
        for i in range(B.shape[0]):
            f.write(str(int(B[i])))
            f.write('\n')
        f.write('\n')
    if (np.array(boundary_array.shape)).shape[0] != 1:
        for i in range(boundary_array.shape[0]):
            a = boundary_array[i][0]
            b = boundary_array[i][1]
            c = boundary_array[i][2]
            B = A[a,b,c,:][A_bool[a,b,c,:]]
            if i == 0:
                C = B
            else:
                C = np.append(C,B)
        C = np.sort(np.unique(C))
        for i in range(C.shape[0]):
            f.write(str(int(C[i])))
            f.write('\n')
        f.write('\n')
    boundary_number += 1
    return boundary_number