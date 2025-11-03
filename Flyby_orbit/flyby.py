import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import BlackHolePO as bh 
import numpy as np

M = 1 
a = np.sqrt(0.84) 
Q = 0.2

E = np.sqrt(2) 
Lz = 4.2
m = 1
q = 0

#intervalo para plot sugerido por E_x_r
N = 10
interval = (E-0.9685)/N
Wlines = []
sp = bh.Space_Time()
names = []
for i in range(N+1):

    f = bh.Linha_de_Mundo([0,13.0, np.pi/2, 0.0,-1,1],[E,Lz,m,q,M,Q,a],Qcarter=0.0, N = 50000,dt = 0.01) 
    sp.append(f)
    names.append("Trajetória :"+str(E))
    E = E - interval
sp.generate_space_time()
sp.plot_WLines(title = "Órbitas não ligadas",names = names , M = M, a = a , Q = Q) 
