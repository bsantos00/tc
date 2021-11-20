# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plot


#Dados ºC, mm
Ts = 80
Tinf = 25
t = 0.5
L = 40
n = 8 #número de nós
l = 2 #número de linhas de nós
#Dados escolhidos

#Material: prata
k = 426 #W/(mm*K)

#convecção forçada em ar
#assume-se:
h = 200 #W/(m^2*K)

#deltax  e delta 
dx= L*l/n/1000
dy= t/2/(l-1)/1000

A=np.zeros((n,n))
C=np.zeros(n)



for j in range(1,n+1):
    match j:
        case 1: 
            A[j-1,j-1]=-k*(dx/dy+dy/dx)-h*dx
            A[j-1,j]=k*dy/(2*dx)
            A[j-1,j+3]=k*dx/dy
            C[j-1]=-h*dx*Tinf-k*dy/(2*dx)*Ts
            
        case (2|3):
            A[j-1,j-1]=-k*(dx/dy+dy/dx)-h*dx
            A[j-1,j]=k*dy/(2*dx)
            A[j-1,j-2]=k*dy/(2*dx)
            A[j-1,j+3]=k*dx/dy
            C[j-1]=-h*dx*Tinf   
            
        case 4:
            A[j-1,j-1]=-1/2*(k*(dx/dy+dy/dx)+h*(dx+dy))
            A[j-1,j-2]=k*dy/(2*dx)
            A[j-1,j+3]=k*dx/(2*dy)
            C[j-1]=-h*1/2*(dx+dy)*Tinf
            
        case 5:
            A[j-1,j-1]=-k*(dx/dy+dy/dx)
            A[j-1,j]=k*dy/(2*dx)
            A[j-1,j-5]=k*dx/dy
            C[j-1]=-k*dy/(2*dx)*Ts
            
        case (6|7): 
            A[j-1,j-1]=-k*(dx/dy+dy/dx)
            A[j-1,j-2]=k*dy/(2*dx)
            A[j-1,j]=k*dy/(2*dx)
            A[j-1,j-5]=k*dx/dy
            
        case 8:
            A[j-1,j-1]=-1/2*(k*(dx/dy+dy/dx)+h*dy)
            A[j-1,j-2]=k*dy/(2*dx)
            A[j-1,j-5]=k*dx/(2*dy)
            C[j-1]=-h*(dy/2)*Tinf  
              

X=np.linalg.inv(A).dot(C)
print(X)

    