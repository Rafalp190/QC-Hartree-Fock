#Integral reading functions obtained from http://www.evangelistalab.org/chem-532-spring-2014/

import math
import sys
import numpy as np 


def read_one_electron_integrals(filename):
    """Read the one-electron integrals contained in the file
       'filename' and return a numpy matrix"""
    ints_lines = open(filename,"r").readlines()
    k = len(ints_lines[0].split())
    ints = np.zeros( (k,k) )
    for i,line in enumerate(ints_lines):
        for j,value in enumerate(line.split()):
            ints[i][j] = float(value)
    return ints


def read_two_electron_integrals(filename,nao):
    """Read the two-electron integrals contained in the file
       'filename' and return a 4-dimensional numpy array that stores
       the integrals in chemistry notation as
       ints[m][n][r][s] = (mn|rs).
       This function needs to know the total number of atomic
       orbitals 'nao'."""
    ints_lines = open(filename,"r").readlines()
    ints = np.zeros( (nao,nao,nao,nao) )
    for line in ints_lines:
        line_split = line.split()
        p,q,r,s = [int(t) for t in line_split[0:4]]
        value = float(line_split[-1])
        ints[p][q][r][s] = ints[q][p][r][s] = value
        ints[p][q][s][r] = ints[q][p][s][r] = value
        ints[r][s][p][q] = ints[r][s][q][p] = value
        ints[s][r][p][q] = ints[s][r][q][p] = value
    return ints