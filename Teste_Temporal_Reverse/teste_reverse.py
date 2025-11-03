import BlackHolePO as bh 
import numpy as np 


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
ln1 = bh.Linha_de_Mundo([0.0,r,np.pi/2.0, 0.0, -1.0, 1.0], [E, Lz, m, q, M, Q, a], dt = 0.1, N = 500, Qcarter = 0)

ln2 = bh.Linha_de_Mundo([0.0,r,np.pi/2.0, 0.0, -1.0, 1.0], [E, Lz, m, q, M, Q, a], dt = 0.1, N = 500, Qcarter = 0,reverse = 1)

sp = bh.Space_Time()
sp.append(ln1)
sp.append(ln2)
sp.generate_space_time()
print(ln2.x)
sp.plot_WLines(title = ' Órbitas circulares em Kerr-Newman ')