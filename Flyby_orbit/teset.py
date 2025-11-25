
import BlackHolePO as bh 
import numpy as np

M = 2.0
a = 0
Q = 0
j = 60
i = 60
r0 = 1000
theta0 = np.pi/2


m = 0
q = 0


delta = r0*r0 + a*a + Q*Q - 2*M*r0
sigma = r0*r0 + a*a*np.cos(theta0)*np.cos(theta0)
g_tt = - (delta - a*a*np.sin(theta0)*np.sin(theta0))/sigma
g_tphi = -a*np.sin(theta0)*np.sin(theta0)*(r0*r0 + a*a - delta)/sigma
g_phiphi = ((r0*r0+a*a)*(r0*r0+a*a) - delta*a*a*np.sin(theta0)*np.sin(theta0))*np.sin(theta0)*np.sin(theta0)/sigma
g_rr = sigma/delta
g_thetatheta = sigma
 
termx = 0.1/(np.sqrt(g_phiphi)*r0)
termy = 0.1/(np.sqrt(g_thetatheta)*r0)
r_bh = M + np.sqrt(M*M -a*a -Q*Q)

Lz = i*termx*g_phiphi + g_tphi
E = g_tphi*termx*i + g_tt
ptheta = j*0.1*np.sqrt(g_thetatheta)/r0
pr = np.sqrt((-g_tt - g_phiphi*termx*termx*i*i - 2*g_tphi*termx*i - g_thetatheta*termy*termy*j*j)*g_rr)

#TESTE
H = -(((-(r0**2 + a**2)*E +a*Lz)**2)/(delta*sigma)) + (((Lz - a*E*np.sin(theta0)*np.sin(theta0))**2)/(sigma*np.sin(theta0)*np.sin(theta0))) + (delta*pr*pr/sigma) + (ptheta*ptheta/sigma)

print("Hamiltoniano:",H)
print("y:", (ptheta*r0/(np.sqrt(g_thetatheta)*pr)))
print("x:", Lz*r0/(np.sqrt(g_phiphi)*pr))
f = bh.Linha_de_Mundo([0,r0, theta0, 0.0,pr,ptheta],[E,Lz,m,q,M,Q,a], N = 500000,dt = 0.01,reverse=1) 
sp = bh.Space_Time()
sp.append(f)
sp.generate_space_time()
sp.plot_Graph([0])
sp.plot_WLines()
