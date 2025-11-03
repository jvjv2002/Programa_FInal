import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from mpl_toolkits.mplot3d import Axes3D
import time
#Fonte: latin1 (Latex)
#Licença: Glise (Liçenca USP)

#ds^2 = -(1-2GM/c^2R) c^2dt^2 + (1-2GM/Rc^2)^(-1) dr^2 + r^2domega^2 

#Constantes notáveis
rs = 5 #2GM/c^2 (Raio de Schwarzschild)
theta0 = np.pi/2 

#Passo temporal: 
h = 0.3

###################################### Realiza transformação de coordenadas #################################

def coordinates(r,phi):
    x = r*np.cos(phi)
    y = r*np.sin(phi)
    return x,y

###################################### Fim da transformação de coordenadas ###################################
###################################### Símbolos de Christoffel ###############################################################
def christoffel(x): 
    # x: Os Símbolos de Christoffel só dependem da posição, pois são características do espaço tempo
    christoffel = np.zeros((4,4,4))
    #dicionário 
    # t -> 0
    t = x[0]
    # r -> 1
    r = x[1] 
    # theta -> 2
    theta = x[2]
    # phi -> 3
    phi = x[3]
    
    christoffel[0] [1][0] = rs/(2*r*(r-rs))
    christoffel[0] [0][1] = christoffel[0] [1][0]#Simetrização 
    
    #(r,0,0)
    christoffel[1] [0][0] = (rs/2)*(r-rs)/(r**3) # 'Força Gravitacional'
    #(r,r,r)
    christoffel[1] [1][1] = -rs/(2*r*(r-rs))
    #(r,theta,theta)
    christoffel[1] [2][2] = -(r-rs)
    #(r,phi,phi)
    christoffel[1] [3][3] = -(r-rs)*((np.sin(theta))**2)
    
    #(theta,r,theta)
    christoffel[2] [1][2] = 1/r
    christoffel[2] [2][1] = christoffel[2] [1][2] #Simetrização 
    #(theta,phi,phi)
    christoffel[2] [3][3] = -np.sin(theta)*np.cos(theta)
    
    #(phi,r,phi)
    christoffel[3] [1][3] = 1/r 
    christoffel[3] [3][1] = christoffel[3] [1][3]
    
    #(phi,theta,phi)
    christoffel[3] [2][3] = np.cos(theta)/np.sin(theta)
    
    return christoffel
#################################################### Fim de Christoffel ##################################################


######################################### Implementação do Método de Runge-Kutta ###############################

    
def dudt(u,x,christoffel):
    
    #Definição da equação a ser solucionada da geodésica (du/dtau)
    du = np.zeros(4)
    for i in range(4):
        for j in range(4):
            for k in range(4):
                du[k] = du[k] - u[i]*u[j]*christoffel[k] [i][j] 

    return du 

    #Atualiza 4-velocidades e posição do evento (Iteração de Runge-Kutta)
def iterate(u,x,dt):# Dada condições iniciais e dt(A discretização temporal)
    
    #K1,L1
    christoff = christoffel(x) #Christoffel para passo inicial
    k1 = dudt(u,x,christoff)
    l1 = u
    
    #K2,L2
    christoff = christoffel(x + (l1*dt/2) ) #Christoffel do primeiro passo intermediário
    k2 = dudt(u + (k1*dt/2) ,x + (l1*dt/2) , christoff)
    l2 = (u + (k1*dt/2))
        
    #K3,L3 
    christoff = christoffel(x + (l2*dt/2) ) #Christoffel do segundo passo intermediário 
    k3 = dudt(u + (k2*dt/2), x + (l2*dt/2) , christoff)  
    l3 = (u + (k2*dt/2))
    
    #K4,L4
    christoff = christoffel(x + l3*dt)
    k4 = dudt(u + k3*dt, x + l3*dt, christoff)
    l4 = (u + k3*dt)
    
    #Atualização das 4-velocidades e posições
    xi = x + dt*(l1 + 2*l2+2*l3 + l4)/6 
    ui = u + dt*(k1 + 2*k2+2*k3 + k4)/6 
    
    return ui,xi  # Retorna 4-velocidades e novas posições
    
    #Implementação do loop:
    
    
    #Parâmetros
dt0 = 0.01
N = 10000
t = [i*dt0 for i in range(N)]
xt = [[0,0,0] for i in range(N)] #(r,theta,phi) para todos os valores possíveis

#Para plot
xc = np.zeros(N) 
yc = np.zeros(N)

#Condições iniciais
r0 = 10
u = np.array([0,0,0, 0.1]) 
x = np.array([0, r0 , theta0 , 0])
u[0] = np.sqrt( (1 + (u[3]**2)*(x[1]**2))/(1-(rs/r0)) )
print("u0:",u[0])
#Energia
E = (1-(rs/r0))*u[0] 
print("Energia:",E)
lz = (r0*r0)*u[3]
print("Lz:",lz)
xt[0] = [x[1],x[2],x[3]]

# Transformação para plot
coords = coordinates(x[1],x[3])
xc[0] = coords[0]
yc[0] = coords[1]
# Loop 

print("###############  Inicializando Força Bruta:  ###############")
print("")
startT = time.time()
for i in range(N-1):
    u,x = iterate(u,x,dt0)
    xt[i+1] = [x[1],x[2],x[3]]
    # Transformação para plot
    coords = coordinates(x[1],x[3])
    xc[i+1] = coords[0]
    yc[i+1] = coords[1]
    
print("Tempo de execução: %s seconds " % (time.time() - startT))
print("")
print("####################################")
print("")
print("")

####################################### Fim de Implementação de Runge Kutta ##################################

####################################### Implementação de Runge-Kutta com Simetria ############################

#   Condições Iniciais
#r0 = 10
#u = np.array([0,0,0, 0.1]) 
#x = np.array([0, r0 , theta0 , 0])
#u[0] = np.sqrt( (1 + (u[3]**2)*(x[1]**2))/(1-(rs/r0)) )

# Simetrias da métrica de Schwarzschild
#e = (1-(rs/x[1]))*u[0]
#l = u[3]*(x[1]*x[1])
#Normalização da quadrivelocidade: uu = -1


#def drdt(r):
    #drdt = sqrt( (-1 - (r^2)uphi^2)  + e^2   )
    #u0 = e/(1-2GM/r)
#    ur = np.sqrt( ((e*e)-1) + rs/r - l*l/(r*r) + rs*l*l/(r*r*r)  )
    
#    return ur

#def iterates(u,x,dt):# Dada condições iniciais e dt (A discretização temporal)
#    xi = np.zeros(4)
#    ui = np.zeros(4)
    
    #Estimativa Inicial
#    k1 = [u[0],drdt(x[1]),u[2],u[3]]
    
    #   Primeiro passo intermediário
#    r2 = x[1] + (k1[1]*dt/2) 
#    k2 = [e/(1-rs/r2) , drdt(r2) , u[2], l/(r2*r2) ] 
    
    #   Segundo passo intermediário
#    r3 = x[1] + (k2[1]*dt/2)
#    k3 = [e/(1-rs/r3) , drdt(r3) , u[2], l/(r3*r3) ]
    
    #   Passo Final
#    r4 = x[1] + (k3[1]*dt)
#    k4 = [e/(1-rs/r4) , drdt(r4) , u[2], l/(r4*r4) ]
    
    #   Finalização
#    xi = x + np.array((k1+ 2*k2 + 2*k3 + k4)/6) #Estimativa Final área dr/dt
#    rf = xi[0]
#    ui = np.array([e/(1-rs/rf) , drdt(rf) , u[2], l/(rf*rf) ])
    
#    return ui,xi

# Implementação do loop
#ts = [i*dt0 for i in range(N)]
#xts = [[0,0,0] for i in range(N)] #(r,theta,phi) para todos os valores possíveis
#xts[0] = [x[1],x[2],x[3]]

#Para plot
#xcs = np.zeros(N) 
#ycs = np.zeros(N)
# Transformação para plot
#coords = coordinates(x[1],x[3])
#xcs[0] = coords[0]
#ycs[0] = coords[1]

# Loop

#print("###############  Inicializando Força Simétrica:  ###############")
#print("")
#startT = time.time()
#for i in range(N-1):
#    u,x = iterates(u,x,dt0)
#    xts[i+1] = [x[1],x[2],x[3]]
#    # Transformação para plot
#    coords = coordinates(x[1],x[3])
#    xcs[i+1] = coords[0]
#    ycs[i+1] = coords[1]
    
#print("Tempo de execução: %s seconds " % (time.time() - startT))
#print("")
#print("####################################")
#print("")
#print("")


#Desenhar espaço-tempo de Schwarzschild
def Sch(x,y): 
    
    r = np.sqrt((x**2 + y**2)) 
    z = 2*np.sqrt(rs*(r-rs)) - 20*rs  # Valor estimado ao se projetar em uma superfície o buraco negro
    return z 

plt.style.use('dark_background')
# Definindo o intervalo dos plots tal qual a sua resolução
x = np.linspace(-20, 20, 50) 
y = np.linspace(-20, 20, 50)

# Criando uma grade de coordenadas
X, Y = np.meshgrid(x, y)

# Definindo a função z = f(x, y)
Z = Sch(X,Y)

# Criando a figura e o eixo 3D
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# Plotando a superfície

surf = ax.plot_surface(X, Y, Z, edgecolor='white',color = 'black',alpha = 0.5,zorder = 1,linewidth = 0.1)

ax.set_title('Buraco Negro de Schwarzschild') # Adicionar título ao gráfico
line, = ax.plot([], [] , [] , color='yellow', label='Trajetória', linewidth=1.0,alpha=1.0,zorder = 2)
#line2, = ax.plot([], [] , [] , color='green', label='Trajetória', linewidth=1.0,alpha=1.0,zorder = 3)
# Remover todos os componentes dos eixos
ax.set_axis_off()
#def update(i):
#    line[1].set_data
#Criar animação 
def init(): #Função Iniciar
    line.set_data([], [])
    line.set_3d_properties([])
#    line2.set_data([], [])
#    line2.set_3d_properties([])
    return line, 

def animate(i): #Atualização frame/frame 
    line.set_xdata(xc[:100*i])
    line.set_ydata(yc[:100*i])
    line.set_3d_properties(Sch(xc[:100*i], yc[:100*i]))
#   line.set_xdata(xcs[:100*i])
#   line.set_ydata(ycs[:100*i])
#   line.set_3d_properties(Sch(xcs[:100*i], ycs[:100*i]))
    return line,
    

ani = animation.FuncAnimation(fig, animate, frames=N//100, init_func=init, interval=100, blit=False)

plt.show() #Mostrar esquema do buraco negro


#Resolver equação diferencial

#Tipo-Luz 


    








