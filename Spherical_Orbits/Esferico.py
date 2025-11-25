import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from scipy.optimize import minimize
import BlackHolePO as bh 
import numpy as np


#Espaço-tempo 
M = 1 
a = 0.91651
Q = 0.25

#Condições Iniciais
m = 1
q = 0

Lz0 = 0 # Imposição do livro -> Black Hole: A student text, third edition
Lz = Lz0

r0 = 5*M

#Black Hole: A student text, third edition
def delta(r): return r**2 + a**2 + Q**2 - 2*M*r
def E_sph(r): 
    E = np.sqrt( r*delta(r)*delta(r)/((r**2 + a**2)*(r**3 - 3*M*r*r + a*a*r + 2*Q*Q*r + a*a*M)) )
    return E

E0 = E_sph(r0)
E = E0
Qcarter = ((E_sph(r0)**2)* (((r0**2+a**2)**2) ) /(delta(r0))) - r0**2 - a*a*E_sph(r0)*E_sph(r0)

# Quantidade de plots
N = 10
sp = bh.Space_Time()
names = []
interval = 1/N
i = 0
for i in range(N+1):
    p1 = bh.Linha_de_Mundo([0,r0,np.pi/2,0,1,1], [E,Lz,m,q,M,Q,a] ,N = 50000 ,dt =0.01,Qcarter = Qcarter)
    names.append("q:"+f'{q:.2f}')
    sp.append(p1)
    q = q + interval
    #   Atualiza energia e momento angular para manter as mesmas condições iniciais mecânica
    #   Isso vale em condições inciais em que \theta = pi/2
    E = E0 + q*Q*r0/(r0**2)
    Lz = Lz0 + q*Q*a*r0/(r0**2)

print("Energia:",E)
print("Constante de Carter",Qcarter)

#Simular equações 

sp.generate_space_time()
index = [i+1 for i in range(N)]
sp.plot_Graph(index)
sp.plot_WLines(title= "Caso 4: Perturbação de movimento esférico", names= names, M = M, a = a, Q = Q, inside= 3)