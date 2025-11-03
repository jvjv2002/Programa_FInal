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
a = 0.5   
Q = 0 
rs = 2*M 
r0 = 4.5*rs
Lz0 = np.sqrt(0.5*rs*r0)*(r0**2 - 2*a*np.sqrt(0.5*rs*r0) + a**2)/((r0*np.sqrt(r0**2 - 1.5*rs*r0 + 2*a*np.sqrt(0.5*rs*r0))))
print("Lz:",Lz0)
#Dados das partículas
marr = np.array([1.0 for i in range(N)])
qarr = np.array([0.0 for i in range(N)])
Lzarr = np.array([Lz0 for i in range(N)])

# Caso de Schwarzschild 
Marr = np.array([ 1.0 for i in range(N)])  
aarr = np.array([ 0.0 for i in range(N)])   
Qarr = np.array([ 0.0 for i in range(N)])  
x_valuesS = np.linspace(0, 100*rs , N)  # Valores de r
y_valuesS = Veff(Marr,Qarr,aarr,marr,qarr, Lzarr, x_valuesS)
# Cria gráfico de Veff

# Caso de Kerr 
Marr = np.array([ 1.0 for i in range(N)])  
aarr = np.array([ 0.5 for i in range(N)])   
Qarr = np.array([ 0.0 for i in range(N)])
x_valuesK = np.linspace(0, 100*rs , N)  # Valores de r
y_valuesK = Veff(Marr,Qarr,aarr,marr,qarr, Lzarr, x_valuesK)
# Cria gráfico de Veff

#Caso de Reissner–Nordström
Marr = np.array([ 1.0 for i in range(N)])  
aarr = np.array([ 0.0 for i in range(N)])   
Qarr = np.array([ 0.5 for i in range(N)])
x_valuesRN = np.linspace(0, 100*rs , N)  # Valores de r
y_valuesRN = Veff(Marr,Qarr,aarr,marr,qarr, Lzarr, x_valuesRN) 
raio = M + np.sqrt(M**2 - a**2 - (0.5)**2)


#Caso de Kerr-Newman 
Marr = np.array([ 1.0 for i in range(N)])  
aarr = np.array([ 0.5 for i in range(N)])   
Qarr = np.array([ 0.5 for i in range(N)])
x_valuesKN = np.linspace(0, 100*rs , N)  # Valores de r
y_valuesKN = Veff(Marr,Qarr,aarr,marr,qarr, Lzarr, x_valuesKN) 
raio = M + np.sqrt(M**2 - a**2 - (0.5)**2)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(x_valuesS,y_valuesS, lw = 0.5)
ax.plot(x_valuesK,y_valuesK, lw = 0.5)
ax.plot(x_valuesRN,y_valuesRN, lw = 0.5)
ax.plot(x_valuesKN,y_valuesKN, lw = 0.5)

plt.legend(["Schwarzschild( M = 1, a = 0, Q = 0)","Kerr(M = 1,a = 0.5, Q = 0 )","Reissner–Nordström(M = 1,a = 0, Q = 0.5 )","Kerr-Newman(M = 1, a = 0.5, Q = 0.5)"], loc='upper right')
ax.set_xlim(raio,30*rs)
ax.set_ylim(0.9, m)


plt.title("Potencial efetivo Veff em função do raio")
plt.show()
 
