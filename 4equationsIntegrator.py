import BlackHolePO as bh 
import numpy as np

# Arquivo para testar possíveis erros numéricos nos polos

#É necessário integrar

# Integra a partir das 4 EDO's de carter que não permite mudança de sinal de r

def dudt(x,const):
    t,r,theta,phi = x 
    m,q,E,Lz,Qcarter,M,a,Q = const 
    P = E*(r**2 + a**2) - a*Lz - q*Q*r 
    delta = r**2 + a**2 + Q**2 - 2*M*r
    sigma = r**2 + (a**2)*((np.cos(theta))**2) 
    print("ERRO?",(-Qcarter - (Lz-a*E)**2 - m*m*r*r + (P*P/delta) )/delta)
    print("ERRO:", P**2 - delta*(m*m*r*r + Qcarter + (Lz - a*E)**2 ))
    dr = (1/sigma)*np.sqrt(P**2 - delta*(m*m*r*r + Qcarter + (Lz - a*E)**2 ))
    dtheta = (1/sigma)*np.sqrt(Qcarter - (np.cos(theta)*np.cos(theta)*(a*a*(m*m-E*E) + ((Lz)*(Lz)/(np.sin(theta)**2)) )) )

    dphi = (1/sigma)*(-(a*E - (Lz/((np.sin(theta))**2)) ) + (a*P/delta))
    dt = (1/sigma)*(-a*(a*E*np.sin(theta)*np.sin(theta) - Lz)   +  ((r*r + a*a)*(P)/delta)    ) 
    
    if(10**3 < (Lz*Lz/((np.sin(theta))**2))):
        print("Possível Erro")
    dudt = np.array([dt,dr,dtheta,dphi])
    
    return dudt
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
# Quantidade de plots
N = 3
sp = bh.Space_Time()
names = []
interval = 1/N
i = 0

#p1 = bh.Linha_de_Mundo([0,r0,np.pi/2,0,1,1], [E,Lz,m,q,M,Q,a] ,Qcarter = Qcarter)
q = q + interval
E = E0 + q*Q*r0/(r0**2)
Lz = Lz0 + q*Q*a*r0/(r0**2)
ln2 = bh.Linha_de_Mundo([0.0, r0, np.pi/2, 0.0], [m,q,E,Lz,Qcarter,M,a,Q],dudt = dudt )

#Testa instabilidad numérica na região
p1 = bh.Linha_de_Mundo([0,r0,np.pi/2,0,1,1], [E,Lz,m,q,M,Q,a] ,Qcarter = Qcarter)
p2  = bh.Linha_de_Mundo([0,r0,np.pi/2,0,1,1], [E,Lz,m,q,M,Q,a] , N = 100000, dt = 0.002,Qcarter = Qcarter)
p3 = bh.Linha_de_Mundo([0,r0,np.pi/2,0,1,1], [E,Lz,m,q,M,Q,a], N = 500000, dt = 0.001 ,Qcarter = Qcarter)
p4 = bh.Linha_de_Mundo([0,r0,np.pi/2,0,1,1], [E,Lz,m,q,M,Q,a], N = 1000000, dt = 0.0005 ,Qcarter = Qcarter)
#names.append("q2:"+f'{q:.2f}')
#sp.append(ln2)
sp.append(p1)
sp.append(p2)
sp.append(p3)
sp.append(p4)
names.append("traj1")
names.append("traj2")
names.append("traj3")
names.append("traj4")
print("Energia:",E)
print("Constante de Carter",Qcarter)

#Simular equações 
#dudt([0.0, r0, np.pi/2, 0.0],[m,q,E,Lz,Qcarter,M,a,Q])
sp.generate_space_time()
sp.plot_Graph([0])
sp.plot_WLines(title= "Movimento esférico", names= names, M = M, a = a, Q = Q, inside= 3)