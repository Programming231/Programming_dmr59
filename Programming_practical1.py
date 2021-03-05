#Hueckel Solver
#1
#hard code a specific 3x3 matrix and find eigenvalues
import numpy as np


v=np.ndarray((3,3))
v=np.array([2,2,4,3,3,5,4,4,6]).reshape(3,3)
y, x=np.linalg.eig(v)

#2
#define a matrix that returns a list of eigenvalues
def get_evals(matrix):      #define function 
    lamda, v=np.linalg.eig(matrix)      #linear algorithm analysis of matrix H
    global evals, evecs
    evals=list(lamda)       #set list of eigenvalues as evals
    evecs=list(v)           #set list of eigen vectors as evecs
    

#3
#write an nxn Huckel matrix for poly-ene
def linearpolyene(n):
    Lh=np.zeros((n,n))              #matrix of zeros (n,n)
    alpha=0
    beta=-1
    for i in range(n):              
        for j in range(n):    
            if i==j:
                Lh[i,j]=alpha           #diagonal elements all alpha
            elif (i-j)**2==1:
                Lh[i,j]=beta      
            else:
                Lh[i,j]=0
    answer=input('Do you want to see the matrix? Yes or No?')
    if answer=='yes':
        print('Matrix:')
        print(Lh)
    get_evals(Lh)        #eigenvalues of nxn huckel matrix
    degeneracy(Lh)

#4
#cyclic
def cyclicpolyene(n):
    Ch=np.zeros((n,n))
    alpha=0
    beta=-1
    for i in range(n):
        for j in range(n):    
            if i==j:
                Ch[i,j]=alpha
            elif (i-j)**2==1 or (i-j)**2==(n-1)**2:
                Ch[i,j]=beta
            else:
                Ch[i,j]=0
    answer=input('Do you want to see the matrix? Yes or No?')
    if answer=='yes':
        print('Matrix:')
        print(Ch)
    get_evals(Ch)        #eigenvalues of nxn huckel matrix
    print(evals)
    degeneracy(Ch)
#5
def degeneracy(matrix):
    get_evals(matrix)
    evals.sort()
    for i in range(len(evals)):
        evals[i]=round(float(evals[i]),3)
    print('Eigenvalue analysis:')
    for i in range(len(evals)):
        v=evals[i-1]
        if evals[i]!=v:
            degeneracy=evals.count(evals[i])  
            print('Degeneracy '+str(degeneracy)+' Eigenvalue: '+str(evals[i]))
        else:
            v=evals[i]
        
    
#6
#sp2 hybridized platonic solids: tetrahedron, cube, dodecahedron
global cube
alpha=0
beta=-1
t=np.array([[0,1,1,1],
            [1,0,1,1],
            [1,1,0,1],
            [1,1,1,0]])
tetrahedron=beta*t+alpha*np.identity(4)

c=np.array([[0,1,1,0,1,0,0,0],
    [1,0,0,1,0,1,0,0],
    [1,0,0,1,0,0,1,0],
    [0,1,1,0,0,0,0,1],
    [1,0,0,0,0,1,1,0],
    [0,1,0,0,1,0,0,1],
    [0,0,1,0,1,0,0,1],
    [0,0,0,1,0,1,1,0]])
cube=beta*c+alpha*np.identity(8)

d=np.array([ [ 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0 ], 
  [ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0 ], 
  [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0 ], 
  [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0 ], 
  [ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1 ], 
  [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0 ], 
  [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0 ], 
  [ 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ], 
  [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ], 
  [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 ], 
  [ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0 ], 
  [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0 ] ])
dodecahedron=beta*d+alpha*np.identity(20)

o=np.array([ [ 0, 1, 1, 1, 1, 0 ], [ 1, 0, 1, 1, 0, 1 ], [ 1, 1, 0, 0, 1, 1 ], 
  [ 1, 1, 0, 0, 1, 1 ], [ 1, 0, 1, 1, 0, 1 ], [ 0, 1, 1, 1, 1, 0 ] ])
octahedron=beta*o+alpha*np.identity(6)

i=np.array([ [ 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0 ], 
  [ 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0 ], 
  [ 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1 ], 
  [ 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1 ], 
  [ 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0 ], 
  [ 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1 ], 
  [ 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0 ], 
  [ 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0 ], 
  [ 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1 ], 
  [ 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 1 ], 
  [ 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0 ], 
  [ 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0 ] ])
icosahedron=beta*i+alpha*np.identity(12)


def platonicsolid(matrix):
    answer=input('Do you want to see the matrix, yes or no?').lower()
    if answer=='yes':
        print('Matrix:')
        print(matrix)
    degeneracy(matrix)


#7
#Into one program
shape=input('Which molecule do you want the Huckel analysis for?').lower()
if shape=='cyclic polyene':
    n=input('How many atoms are in this cyclic polyene?')
    cyclicpolyene(int(n))
if shape=='linear polyene':
    n=input('How many atoms are in this linear polyene?')
    linearpolyene(int(n))
if shape=='cube':
    platonicsolid(cube)
if shape=='tetrahedron':
    platonicsolid(tetrahedron)
if shape=='dodecahedron':
    platonicsolid(dodecahedron)
if shape=='icosahedron':
    platonicsolid(icosahedron)
if shape=='octahedron':
    platonicsolid(octahedron)
