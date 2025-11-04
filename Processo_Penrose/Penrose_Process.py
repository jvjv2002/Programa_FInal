import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import BlackHolePO as bh 
import numpy as np
M = 3
a = 2.98
Q = 0

#Ponto de virada para 0 
r = 1.80 * M
def r_erg(theta): return (M + np.sqrt(M**2 - (a*np.cos(theta))**2 - Q**2))
r_bh = M + np.sqrt(M**2 - a**2 - Q**2) 
print("Ponto de ocorrencia da separação na ergosfera:",r)
print("r da ergosfera", r_erg(np.pi/2))
print("r do buraco negro", r_bh)
delta = r**2 + a*a + Q*Q - 2*M*r

#Partícula 0: Partícula maciça inicial
E0 = 1
Lz0 = (-2*a*M + np.sqrt(2*M*r*delta))/(r-2*M)
m0 = 1 
Q0 = 0 #No equador

    #Linha de mundo  
p0 = bh.Linha_de_Mundo([0,r,np.pi/2,0,-1,1],[E0,Lz0,m0,0,M,Q,a],Qcarter = 0, reverse = 1) #Queremos ver por onde ela passou 
print("Condições Inciais",p0.x)


#Partícula 1: Fóton que cai no buraco negro
E1 = -0.5*(np.sqrt(2*M/r) -1 )
Lz1 = (-2*a*M - r*np.sqrt(delta))*E1/(r-2*M)
m1 = 0 # Se trata de um fóton 
Q1 = 0 # No equador 
    #Linha de mundo 
p1 = bh.Linha_de_Mundo([0,r,np.pi/2,0,-1,1],[E1,Lz1,m1,0,M,Q,a],N = 2150,Qcarter = 0)
print("Condições Inciais:",p1.x)

#Partícula 2: Fóton que escapa do buraco negro
E2  = 0.5*(np.sqrt(2*M/r) + 1)
Lz2 = (-2*a*M + r*np.sqrt(delta))*E2/(r-2*M)
m2 = 0 # Se trata de um fóton 
Q2 = 0 # No equador

    #Linha de mundo
p2 = bh.Linha_de_Mundo([0,r,np.pi/2,0,-1,1],[E2,Lz2,m2,0,M,Q,a],Qcarter = 0)
print("Condições Iniciais:",p2.x)

# Hora de ilustrar o processo 
sp = bh.Space_Time() 
sp.append(p0)
sp.append(p1)
sp.append(p2)
sp.generate_space_time()
index = [0,1,2]
sp.plot_Graph(index)
sp.plot_WLines(title = "Caso de Estudo 4:Processo de Penrose ",names = ["Partícula 0: E="+f'{E0:.2f}'+"; Lz ="+f'{Lz0:.2f}', "Fóton 1: E="+f'{E1:.2f}'+"; Lz ="+f'{Lz1:.2f}', "Fóton 2: E="+f'{E2:.2f}'+"; Lz ="+f'{Lz2:.2f}'],M = M, Q = Q, a = a, inside = 3)

