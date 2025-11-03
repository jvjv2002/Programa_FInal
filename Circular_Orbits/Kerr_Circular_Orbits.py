import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
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
Q = 0.0 
rs = 2*M 
r = 4*rs 
m = 1.0 
q = 0.0

E = (r**2 - rs*r + a*np.sqrt(0.5*rs*r))/(r*np.sqrt(r**2 - 1.5*rs*r + 2*a*np.sqrt(0.5*rs*r)))
Lz = np.sqrt(0.5*rs*r)*(r**2 - 2*a*np.sqrt(0.5*rs*r) + a**2)/((r*np.sqrt(r**2 - 1.5*rs*r + 2*a*np.sqrt(0.5*rs*r))))
print("Lz:",Lz)
Q = 0
print("r+",M + np.sqrt(M**2 - a**2 ))
ln1 = bh.Linha_de_Mundo([0.0,r,np.pi/2.0, 0.0, -1.0, 1.0], [E, Lz, m, q, M, Q, a], dt = 0.1, N = 20000, Qcarter = 0)
ln2 = bh.Linha_de_Mundo([0.0,r,np.pi/2.0, 0.0, -1.0, 1.0], [(1+0.001)*E, (1+0.001)*Lz, m, q, M, Q, a], dt = 0.1, N = 20000, Qcarter = 0)
ln3 = bh.Linha_de_Mundo([0.0,r,np.pi/2.0, 0.0, -1.0, 1.0], [(1+0.01)*E, (1+0.01)*Lz, m, q, M, Q, a], dt= 0.1, N = 20000, Qcarter = 0)
sp = bh.Space_Time()
sp.append(ln1)
sp.append(ln2)
sp.append(ln3)
sp.generate_space_time()
titulo1 = "E = "+f"{E:05f}"+"; Lz = "+f"{Lz:05f}"
titulo2 = "E = "+f"{(1+0.001)*E:05f}"+"; Lz = "+f"{(1+0.001)*Lz:05f}"
titulo3 = "E = "+f"{(1+0.01)*E:05f}"+"; Lz = "+f"{(1+0.01)*Lz:05f}"
index = [i for i in range(3)]
sp.plot_Graph(index)
sp.plot_WLines(title = ' Órbitas circulares em Kerr-Newman ',names = [titulo1,titulo2,titulo3],M = M, a = a, Q = Q)
