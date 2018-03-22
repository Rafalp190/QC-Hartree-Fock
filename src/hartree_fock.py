from integral_reader import *
from numpy import linalg as LA
import numpy as np

def form_G(D, V, nao):
	G = np.zeros((nao, nao))
	for u in range(0,nao):
		for v in range(0,nao):
			for p in range(0,nao):
				for o in range(0,nao):
					G[u][v] = D[p][o] * ((2*V[u][v][p][o]) - V[u][p][v][o])
	return G

def form_E(D, H, F, nao):
	E_temp = 0
	for m in range(0,nao):
		for n in  range(0,nao):
			E_temp += D[m][n]*(H[m][n] + F[m][n])
	return E_temp

def form_D(C, Nhalf, nao):
	D_temp = np.zeros((nao, nao))
	for m in range(0,nao):
		for n in range(0,nao):
			for i in range(0,Nhalf):
				D_temp[m][n] = C[m][i] * C[n][i]
	return D_temp

def form_delta_D(D, D_i, nao):
	delta_D = 0 
	for m in range(0,nao):
		for n in range(0,nao):
			delta_D += (D_i[m][n] - D[m][n])**2
	return delta_D**1/2
# Define the number of atomic orbitals (nao)
nao = 7

# Define the number of electrons
N = 10

# Define the number of filled spatial orbitals
Nhalf = int(N / 2)

VNN = float(open("../static/vnn","r").readlines()[0])
S = read_one_electron_integrals("../static/overlap")
H = read_one_electron_integrals("../static/one-electron")
V = read_two_electron_integrals("../static/two-electron",nao)

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
#print(L)
s_diag = np.diag(s**(-0.5))
X = L @ s_diag @ np.transpose(L)

print("X matrix\n",X)
#Y = np.transpose(X) @ S @ X
#print(Y)


#Initialization of values
D = np.zeros((nao, nao))

k = 0
k_bound = 1
E = 0
non_convergence = True
print("------------------------------------")

while k<k_bound or non_convergence == False:
	G = form_G(D, V, nao)

	print("G matrix")
	print(G)
	#print("H")
	#print(H)

	F = np.add(H, G)
	print("F matrix")
	print(F)
	E_i = VNN + form_E(D, H, F, nao)
	#print(E_i)

	F_i = X @ F @ X
	print("F' matrix")
	print(F_i)
	F_eig = LA.eig(F_i)
	print("Eigenvalues of F'")
	print(F_eig[0])
	print("C'")
	print(F_eig[1])
	#12
	C = X @ F_eig[1]
	#print(C)
	#13
	D_i = form_D(C, Nhalf, nao)
	#print(D_i)
	#14
	
	if k == 0:
		E =  E_i
		#print(E)
		D = D_i
		#print(D)
		k += 1
	else:
		delta_E = E_i - E
		#print(delta_E)
		delta_D = form_delta_D(D, D_i, nao)
		if abs(delta_E) < 10e-9 and abs(delta_D) < 10e-5:
			non_convergence = False
		E =  E_i

		D = D_i
		k += 1
print(k)
print(E)
print(D)
