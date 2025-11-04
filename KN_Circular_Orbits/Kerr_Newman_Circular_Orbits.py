import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scipy.optimize import minimize
import BlackHolePO as bh 
import numpy as np
#ln1 = bh.Linha_de_Mundo( [0.0,10.0, 1.57 ,0.0 ,0.0 ,0.0],  [2.0,20.0, 2 , 0.0 , 2.5 ,0 ,0])
#ln2 = bh.Linha_de_Mundo( [0.0,10.0, 1.57 ,0.0 ,0.0 ,0.0],  [1.4142135623730951 ,10.0 , 1.#0 , 0.0 , 0.0 ,0 ,0])
#sp = bh.Space_Time()
#sp.append(ln1)
#sp.append(ln2)
#sp.generate_space_time() 
#sp.plot_WLines(title='teste')    
M = 1.0 
a = 0.5 
Q = 0.5
rs = 2*M 
r0 = 4.5*rs
m = 1.0 
q = 0.0

Lz = 3.45
print(Lz)
#Cria espaço-tempo
sp = bh.Space_Time()

def Veff(r):
    delta = r**2 + Q**2 + a**2 - 2*M*r
    Sigma = r**2 # Em theta = pi/2
    A = (r**2 + a**2)**2 - delta*a*a
    B = (r**2 + a**2)*(a*Lz + Q*q*r) - a*delta*Lz 
    C = (a*Lz + q*Q*r)**2 - delta*Lz*Lz - delta*Sigma*m*m
    Veff = (B + np.sqrt(B**2 - A*C))/(A) 
        
    return Veff

# Procura Trajetória circular
result = minimize(Veff, 8.4 , tol= 10**(-30) ) 
E0 = result.fun
r0 = result.x[0] 
Lz0 = Lz
E = E0

# Verifica perturbação gerada pela carga da partícula
N = 5
names = []
intervalo = 1/N
for i in range(N+1):
    ln1 = bh.Linha_de_Mundo([0,r0,np.pi/2,0,1,1],[E,Lz,m,q,M,Q,a],N = 500000 ,Qcarter = 0)
    sp.append(ln1)
    names.append("q = "+f'{q:.2f}')
    q = q + intervalo
    E = E0 + q*Q*r0/(r0**2)
    Lz = Lz0 + q*Q*a*r0/(r0**2)
    #Garantir que energia cinéitca+ repouso da partícula se mantêm os mesmos
    #print("Condição Inicial: ",i,ln1.x). De fato, desconderando o campo eletromagnético Condição Inicial mecãnica se manteve
    

sp.generate_space_time()
index = [i for i in range(N)]
sp.plot_Graph(index)

sp.plot_WLines(title = "Caso 2: Órbita circulares perturbada por cargas",names = names ,M = M , Q = Q , a = a, inside = 2)
