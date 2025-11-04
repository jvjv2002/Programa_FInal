import BlackHolePO as bh 
import numpy as np

# É possível resolver as 4 EDO's que envolve quadrados demonstrado por Carter, desde que se escolha os sinais de $R(r)$ e $\Theta(\theta)$

# Testa se existem diferenças numéricas, o que não foi detectado, portanto se existe algum erro numérico devido aos polos, ele não foi detectado a partir da diferença dos dois algoritimos

#É necessário integrar

# Integra a partir das 4 EDO's de carter que não permite mudança de sinal de r
def dudt(x,const):
    t,r,theta,phi = x 
    m,q,E,Lz,Qcarter,M,a,Q = const 
    P = E*(r**2 + a**2) - a*Lz - q*Q*r 
    delta = r**2 + a**2 + Q**2 - 2*M*r
    sigma = r**2 - (a**2)*((np.cos(theta))**2) 
    
    dr = (1/sigma)*np.sqrt(P**2 - delta*(m*m*r*r + Qcarter + (Lz + a*E)**2 ))
    dtheta = (1/sigma)*np.sqrt(Qcarter - (np.cos(theta)*np.cos(theta)*(a*a*(m*m-E*E) + ((Lz)*(Lz)/(np.sin(theta)**2)) )) )
    
    dphi = (1/sigma)*(-(a*E - (Lz/((np.sin(theta))**2)) ) + (a*P/delta))
    dt = (1/sigma)*(-a*(a*E*(np.sin(theta)**2) - Lz ) + (r**2 + a**2)*P/delta)
    
    #print("dtheta/dt",dtheta/dt)
    
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
for i in range(N+1):
    p1 = bh.Linha_de_Mundo([0,r0,np.pi/2,0,1,1], [E,Lz,m,q,M,Q,a] ,Qcarter = Qcarter)
    ln2 = bh.Linha_de_Mundo([0.0,r0,np.pi/2,0,1,1], [E,Lz,m,q,M,Q,a], Qcarter = Qcarter)
    names.append("q:"+f'{q:.2f}')
    names.append("q2:"+f'{q:.2f}')
    sp.append(p1)
    sp.append(ln2)
    q = q + interval
    #Atualiza energia e momento angular para manter as mesmas condições iniciais mecânica
    E = E0 + q*Q*r0/(r0**2)
    Lz = Lz0 + q*Q*a*r0/(r0**2)

print("Energia:",E)
print("Constante de Carter",Qcarter)

#Simular equações 

sp.generate_space_time()
index = [i+1 for i in range(2*N)]
sp.plot_Graph(index)
sp.plot_WLines(title= "Movimento esférico", names= names, M = M, a = a, Q = Q, inside= 3)
