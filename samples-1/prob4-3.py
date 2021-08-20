import math
import numpy as np

n = 5
loop = 5


a = np.matrix([[0, ]*n, ]*n, dtype=np.int64)
for i in range(n-1):
    a[i, i] = 2
    a[i, i+1] = a[i+1, i] = -1
a[n-1, n-1] = 1
print('Matrix A:')
print(a)

# def mul(u, v):
#     ret = np.copy(u)
#     for i in range(n):
#         for j in range(n):
#             ret[i, j] = 0
#             for k in range(n):
#                 ret[i, j] += u[i, k]*v[k, j]
#     return ret


# a2 = np.copy(a)
# a4 = np.copy(a)
# for l in range(loop):
#     for i in range(n):
#         for j in range(n):
#             a4[i, j] = 0
#             for k in range(n):
#                 a4[i, j] += a2[i, k]*a2[k, j]
#     a2 = np.copy(a4)
# print(a2)
# print(a**32)
a2 = np.copy(a**int(2**loop))
v = a2*np.matrix([[1, ]]*n, dtype=np.int64)
print(v)
v2 = np.copy(a*v)
print(v2)
# print(v2.T)
# print(v2.T)
# m = np.copy(v2.T@v2)
# d = np.copy(v2.T@v)
# print(v2.T)
# print(m, d)
# print(m[0, 0]/d[0, 0])
m = 0
d = 0
for i in range(n):
    t = int(v[i])
    t2 = int(v2[i])
    m += 0+t2*t2
    d += 0+t*t2
print(m, d, m/d)
print(2*(1-math.cos((2*n-1)*math.pi/(2*n+1))))

def eigv(k):
  return 2*(1-math.cos((2*k-1)*math.pi/(2*n+1)))

print(1/math.log(eigv(5)/eigv(4)))