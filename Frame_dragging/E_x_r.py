import numpy as np 
import matplotlib.pyplot as plt


N = 1000 # Numero de pontos a serem plotados

def Veff(M,Q,a,m,q,Lz,r):
    delta = r**2 + Q**2 + a**2 - 2*M*r
    Sigma = r**2 # Em theta = pi/2
    A = (r**2 + a**2)**2 - delta*a*a
    B = (r**2 + a**2)*(a*Lz + Q*q*r) - a*delta*Lz 
    C = (a*Lz + q*Q*r)**2 - delta*Lz*Lz - delta*Sigma*m*m
    Veff = (B + np.sqrt(B**2 - A*C))/(A) 
    print(Veff)
        
    return Veff




# E e Lz para quais os gráficos serão analisados
m = 1.0
M = 1.0
a = 0.84
Q = 0.2
rs = 2*M 
q = 0.1
Lz0 = 0.0
Lz1 = -4.04
Lz2 = 4.0
print("Lz:",Lz0)
#Dados das partículas
marr = np.array([1.0 for i in range(N)])
qarr = np.array([q for i in range(N)])
Lzarr = np.array([Lz0 for i in range(N)])
Lzarr1 = np.array([Lz1 for i in range(N)])
Lzarr2 = np.array([Lz2 for i in range(N)])
# Cria gráfico de Veff


Marr = np.array([ M for i in range(N)])  
aarr = np.array([ a for i in range(N)])   
Qarr = np.array([ Q for i in range(N)])
#Partícula 1: Lz = 0
x_valuesK = np.linspace(0, 100*rs , N)  # Valores de r
y_valuesK = Veff(Marr,Qarr,aarr,marr,qarr, Lzarr, x_valuesK)
#Partícula 2: Lz = -4
y_valuesK1 = Veff(Marr,Qarr,aarr,marr,qarr, Lzarr1, x_valuesK)
#Partícula 3: Lz = 4
y_valuesK2 = Veff(Marr,Qarr,aarr,marr,qarr, Lzarr2, x_valuesK)

# Cria gráfico de Veff
fig = plt.figure()
ax = fig.add_subplot(111)
N = 3
colors = plt.cm.gist_rainbow(np.linspace(0, 1, N, endpoint=False))
ax.plot(x_valuesK,y_valuesK, lw = 0.5, color = colors[0])
ax.plot(x_valuesK, y_valuesK1, lw = 0.5, color = colors[1])
ax.plot(x_valuesK,y_valuesK2, lw = 0.5, color = colors[2])
ax.legend(["Partícula 1: Lz = 0","Partícula 2: Lz = -4","Partícula 3: Lz = 4"])
ax.set_xlim(2*rs,30*rs)
ax.set_ylim(0.9, m)
ax.set_xlabel("r/M")
ax.set_ylabel("Veff(r)/M")


plt.title("Potencial efetivo Veff em função do raio")
plt.show()
 
