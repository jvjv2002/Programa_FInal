import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import BlackHolePO as bh 
import numpy as np 


M = 1.0 
a = 0.99
Q = 0
m = 1.0
E = 0.9733232697
Lz = 2.3608795160
r0 = 1.8421052632
zoomWhirl = bh.Linha_de_Mundo([0,r0,np.pi/2,0.0,1.0,1.0],[E,Lz,m,0,M,Q,a], N = 200000, Qcarter = 0)
sp = bh.Space_Time()
sp.append(zoomWhirl)
sp.generate_space_time()
sp.plot_Graph([0])
sp.plot_WLines(M = M ,a = a, Q = Q)