# グラフ
import numpy as np
import matplotlib.pyplot as plt

# Define control points:CP
CP = np.array([[1, 1],
               [1, 2],
               [2, 3],
               [3, 2],
               [5, 2],
               [6, 3],
               [5, 4]])

# Define polynomial order:n
n = 2   # n次のB-スプライン曲線

l = CP.shape[0]  # 制御点の個数:l

m = l + n + 1   # 制御点の個数(l) = m - n - 1 ，m:ノットの個数

# Define knot vector
knot = np.zeros((m))
for i in range(m - 2 * (n + 1)):
    knot[i+n+1] = (i + 1) / (m - 2 * (n + 1) + 1)
for i in range(n + 1):
    knot[m-i-1] = 1.0
print("knot vector =")
print(knot)

# C0連続(n=2)
# knot = np.array([0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0])
# print(knot)

# 描画
fig, ax = plt.subplots()

# Define color vector
color = np.array(["r", "g", "b", "c", "m", "y", "k"])

N = np.zeros((l, n+1))   # スプライン基底関数(l*n+1)
delta = 300             # nowknotの刻み幅
for i in range(delta):
    nowknot = (i/delta) * knot[m-1]  # 1 / 1000ずつ増加
    for k in range(n+1):
        for j in range(l):
            if k == 0:
                if (knot[j] <= nowknot) and (nowknot < knot[j+1]):
                    N[j][k] = 1.0
                else:
                    N[j][k] = 0.0
            else:
                if (knot[j+k] - knot[j]) == 0:
                    a = 0
                if (knot[j+k] - knot[j]) != 0:
                    a = ((nowknot - knot[j]) /
                         (knot[j+k] - knot[j])) * N[j][k-1]
                if (knot[j+k+1] - knot[j+1]) == 0:
                    b = 0
                if (knot[j+k+1] - knot[j+1]) != 0:
                    b = ((knot[j+k+1] - nowknot) /
                         (knot[j+k+1] - knot[j+1])) * N[j+1][k-1]
                N[j][k] = a + b
    Cx = 0
    Cy = 0
    for j in range(l):
        Cx += N[j][n] * CP[j][0]
        Cy += N[j][n] * CP[j][1]
    ax.scatter(Cx, Cy, c=color[4], s=0.1)
for i in range(l):
    x = CP[i][0]
    y = CP[i][1]
    ax.scatter(x, y, c=color[2], s=5)
plt.show()
