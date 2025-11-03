import numpy as np 
import matplotlib.pyplot as plt


Np = 1000 # Numero de pontos a serem plotados

m = 1.0
q = 0.0
M = 1.0
a = 0.5   
Q = 0.5

rs = 2*M 
r0 = 4.5*rs
# E e Lz para quais os gráficos serão analisados
Lz = np.sqrt(0.5*rs*r0)*(r0**2 - 2*a*np.sqrt(0.5*rs*r0) + a**2)/((r0*np.sqrt(r0**2 - 1.5*rs*r0 + 2*a*np.sqrt(0.5*rs*r0))))

print("Lz:",Lz)

def Veff(r):
    delta = r**2 + Q**2 + a**2 - 2*M*r
    Sigma = r**2 # Em theta = pi/2
    A = (r**2 + a**2)**2 - delta*a*a
    B = (r**2 + a**2)*(a*Lz + Q*q*r) - a*delta*Lz 
    C = (a*Lz + q*Q*r)**2 - delta*Lz*Lz - delta*Sigma*m*m
    Veff = (B + np.sqrt(B**2 - A*C))/(A) 
        
    return Veff


N = 5
intervalo = 1/N
fig = plt.figure()
ax = fig.add_subplot(111)
raio = M + np.sqrt(M**2 - a**2 - (Q)**2)
print("r_bh",raio)
names = []
for i in range(N+1):
#Caso de Kerr-Newman 
    Marr = np.array([ M for i in range(Np)])  
    aarr = np.array([ a for i in range(Np)])   
    Qarr = np.array([ Q for i in range(Np)])
    marr = np.array([m for i in range(Np)])
    qarr = np.array([q for i in range(Np)])
    Lzarr = np.array([Lz for i in range(Np)])
    x_valuesKN = np.linspace(0, 100*rs , Np)  # Valores de r
    y_valuesKN = Veff(x_valuesKN) 
    ax.plot(x_valuesKN,y_valuesKN, lw = 0.5)
    names.append("Trajetória: "+str(i)+";  q: "+str(q))
    q = q + intervalo
    
ax.set_xlim(raio,30*rs)
ax.set_ylim(0.9, 1.2*m)
plt.legend(names, loc='upper right')
plt.title("Potencial efetivo Veff em função do raio")
plt.show()
 
