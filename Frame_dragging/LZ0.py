import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import BlackHolePO as bh 
import numpy as np

M = 1.0
a = 0.84
Q = 0.2
r0 = 60.0
r_bh = M + np.sqrt(M**2 - a**2 - Q**2) # r do buraco negro

theta0 = np.pi/2 
def sigma(r,theta): return r**2 + (a*np.cos(theta))**2


#Particula 1: Lz = 0 
m = 1
q = 0.0
E = 0.988 #+ q*Q*r0/sigma(r0,theta0)
Lz = 0 #+ q*Q*a*r0*((np.sin(theta0))**2)/sigma(r0,theta0)
p1 = bh.Linha_de_Mundo([0,r0,np.pi/2,0,-1,1],[E,Lz,m,q,M,Q,a],Qcarter = 0)


#Particula 2 : Lz = -4
m = 1 
q = 0.0
E = 0.988
Lz = -4.0
p2 = bh.Linha_de_Mundo([0,r0,np.pi/2,0,-1,1],[E,Lz,m,q,M,Q,a], N = 100000 ,Qcarter = 0)

#Particula 3 : Lz = 4
m = 1 
q = 0.0
E = 0.988
Lz = 4.0
p3 = bh.Linha_de_Mundo([0,r0,np.pi/2,0,-1,1],[E,Lz,m,q,M,Q,a], N = 100000 ,Qcarter = 0)

sp = bh.Space_Time()
sp.append(p1)
sp.append(p2)
sp.append(p3)
sp.generate_space_time()

sp.plot_Graph([0,1,2])
sp.plot_WLines(title = "Caso de estudo 2: Captura de partículas ",names = ["Partícula 1: Lz = 0","Particula 2: Lz = -4","Particula 3: Lz = 4"], M =M , a=a, Q =Q , inside =1)

ln = sp.returnWline(0)
angular = (ln.xt[-2][3] - ln.xt[-3][3])/(ln.xt[-2][0] - ln.xt[-3][0])
print("----------Toda Partícula Tende a velocidade angular do buraco negro----------")
print("Velocidade Angular do Buraco Negro",a/(r_bh**2 + a**2))
print("Velocidade Angular calculada para p0", angular)
ln = sp.returnWline(1)
angular = (ln.xt[-2][3] - ln.xt[-3][3])/(ln.xt[-2][0] - ln.xt[-3][0])
print("Velocidade Angular calculada para p1", angular)