from integral_reader import *
from numpy import linalg as LA
import numpy as np

# Define the number of atomic orbitals (nao)
nao = 7

# Define the number of electrons
N = 10

# Define the number of filled spatial orbitals
Nhalf = N / 2

VNN = float(open("../lib/vnn","r").readlines()[0])
S = read_one_electron_integrals("../lib/overlap")
H = read_one_electron_integrals("../lib/one-electron")
V = read_two_electron_integrals("../lib/two-electron",nao)

np.set_printoptions(precision=5,linewidth=100)

#print("S matrix:\n",S)
#print("H matrix:\n",H)

# Print the two-electron integrals in chemist notation
# for m in range(0,nao):
#     for n in range(0,nao):
#         for r in range(0,nao):
#             for s in range(0,nao):
#                 print('({:d} {:d}|{:d} {:d}) = {:20.12f}'.format(m, n, r, s, V[m][n][r][s]))

#4 
# a = 0
# b = 0
# for m in range(0,nao):
#     for n in range(0,nao):
#         for r in range(0,nao):
#             for s in range(0,nao):
#             	a = a + (V[m][n][r][s])**2
#             	b = b + abs(V[m][n][r][s])
# print(a)
# print(b)

eigS = LA.eig(S)
s = eigS[0]
L = eigS[1]

s_diag = np.diag(s**(-0.5))
X = L * s_diag * L.transpose()
Y = X.transpose() * S * X

print(Y)