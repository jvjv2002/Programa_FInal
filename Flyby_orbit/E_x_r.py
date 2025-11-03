import numpy as np 
import matplotlib.pyplot as plt


N = 1000 # Numero de pontos a serem plotados

def Veff(M,Q,a,m,q,Lz,r):
    delta = r**2 + Q**2 + a**2 - 2*M*r
    Sigma = r**2 # Em theta = pi/2
    A = (r**2 + a**2)**2 - delta*a*a
    B = (r**2 + a**2)*(a*Lz + Q*q*r) - a*delta*Lz 
    C = (a*Lz + q*Q*r)**2 - delta*Lz*Lz - delta*Sigma*m*m
    Veff = (B + np.sqrt(B**2 - 4*A*C))/(2*A) 
    print(Veff)
        
    return Veff




# E e Lz para quais os gráficos serão analisados
m = 1.0
M = 1.0
a = np.sqrt(0.84)
Q = 0.3
rs = 2*M 
Lz0 = 4.2
print("Lz:",Lz0)
#Dados das partículas
marr = np.array([1.0 for i in range(N)])
qarr = np.array([0.0 for i in range(N)])
Lzarr = np.array([Lz0 for i in range(N)])

# Cria gráfico de Veff

# Caso de Kerr 
Marr = np.array([ M for i in range(N)])  
aarr = np.array([ a for i in range(N)])   
Qarr = np.array([ Q for i in range(N)])
x_valuesK = np.linspace(0, 100*rs , N)  # Valores de r
y_valuesK = Veff(Marr,Qarr,aarr,marr,qarr, Lzarr, x_valuesK)
# Cria gráfico de Veff


fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(x_valuesS,y_valuesS, lw = 0.5)
ax.plot(x_valuesK,y_valuesK, lw = 0.5)
#ax.plot(x_valuesKN,y_valuesKN, lw = 0.5)

ax.set_xlim(2*rs,30*rs)
ax.set_ylim(0.9, m)


print("r_K",min(y_valuesK))

plt.title("Potencial efetivo Veff em função do raio")
plt.show()
 
