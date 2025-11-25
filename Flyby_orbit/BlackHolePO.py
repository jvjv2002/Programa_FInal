import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
from matplotlib.animation import PillowWriter
from mpl_toolkits.mplot3d import Axes3D
from typing import List
import time


class Linha_de_Mundo: # Linha de mundo de uma partícula teste no espaço-tempo de KN
    
    def coordinates(self):
        #Obtêm representação fora de escala em coordenadas cartesianas
        r = self.x[1]
        theta = self.x[2]
        phi = self.x[3]
        
        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)
    
        return(x,y,z)

    def dudt0(x,const):
        # Dinâmica do espaço-tempo de Kerr-Newman 
        
        t, r , theta , phi , pr , ptheta = x
        E, lz,m,q,M,Q,a = const
        # Expressões algébricas úteis
        delta = r**2 + a*a + Q*Q - 2*M*r 
        sigma = r**2 + a*a*(np.cos(theta)**2)
        util = ((a*a + r*r)*E - a*lz - q*Q*r ) # = P
        cosseno = np.cos(theta)
        seno = np.sin(theta)
        
        #dt
        dt = (1/sigma)*(-a*(a*E*seno*seno - lz)   +  ((r*r + a*a)*(util)/delta)    ) 
        
        #dr
        dr = (delta/sigma) * pr 
        #dr = (np.sqrt( ((E*E)-1) -     (lz*lz/(r*r)) + (2*M/r) + (2*M*lz*lz/(r*r*r))          ))

        #dtheta 
        dtheta = (1/sigma) * ptheta
        
        #dphi 
        dphi = (1/sigma)*( -(a*E - lz/(seno*seno)) + (a/delta)*util  )

        #dpr
        dpr = (1/(2*sigma))*(-((util*util*(-2*M+2*r))/(delta*delta)) + (2*util*(2*r*E - q*Q)/(delta)) - ((-2*M + 2*r)*pr*pr) - (2*m*m*r)    )
        
        #dptheta
        dptheta =  (1/sigma)*( cosseno*seno*(a*a*(-E*E + m*m) + lz*lz/(seno*seno)) + ((cosseno**3)*lz*lz/(seno**3)) )
        
        dudt = np.array([dt,dr,dtheta,dphi,dpr,dptheta])

        return dudt #x = ( t, r , theta , phi , pr , utheta)
    
        
    def __init__(self, x ,const , dudt = None, dt= 0.01, N = 50000, Qcarter = None, reverse = 0) :
        """
        Entrada: 
        
        x:float[] -> Condições iniciais de variáveis dinâmicas. Em Kerr-Newman:[t,r,theta,phi,pr,ptheta]
        
        const:float[] -> Constantes de Movimento da partícula teste e parâmetros do espaço-tempo [E,lz,m,q,M,Q,a] OU constantes personalizadas necessárias para dinâmica 
        
        dudt: function -> Dita a dinâmica do sistema 
        
        dt: float -> Tamanho de passo temporal
        
        N: int -> Inteiro para determinar a quantidade de passos temporais para qual o sistema será evoluído
        
        Q: float -> Constante de Carter, se não especificada, considera-se que p_r e p_theta foram dados, caso contrário
        calcula-se eles a partir do mesmo, sendo somente necessário colocar o sinal em x
        
        reverse: int -> 0 se a evolução for para frente no tempo, 1 se a evolução for realizada para o passado
        """
        self.reverse = reverse
        self.x = x # Condições iniciais de variáveis dinâmicas. Em Kerr-Newman:[t,r,theta,phi,pr,ptheta]
        self.const = const # Constantes de Movimento da partícula teste e parâmetros do espaço-tempo [E,lz,m,q,M,Q,a] OU constantes personalizadas necessárias para dinâmica 
        self.xt = [] # História das variáveis dinâmicas que compõem a linha-de-mundo
        self.xplot = []
        self.yplot = []
        self.zplot = []
        # Instruções na discretização da trajetória 
        self.dt = dt 
        self.N = N 
        # Para testes de Captura de partícula-teste
        if(dudt == None):
            M = self.const[4]
            Q = self.const[5]
            a = self.const[6]
            self.r_horizon = M + np.sqrt(M**2 - a**2 - Q**2)
        else:
            self.r_horizon = 0.0
            
            
        #Utiliza a maneira opcional de dar as condições iniciais por meio da Constante de Carter K 
        if(Qcarter != None):
            # Propriedades da partícula
            E = self.const[0]
            lz = self.const[1]
            m = self.const[2]
            q = self.const[3]
            
            # Propriedades do espaço-tempo
            r = self.x[1]
            theta = self.x[2]
            M = self.const[4]
            Q = self.const[5]
            a = self.const[6]
            
            delta = r**2 + a*a + Q*Q - 2*M*r 
            print("delta:",delta)
            sigma = r**2 + a*a*(np.cos(theta)**2)
            util = ((a*a + r*r)*E - a*lz - q*Q*r ) # = P
            print("util",util)
            cosseno = np.cos(theta)
            seno = np.sin(theta)
            
            # A escolha de sinal é necessária em ambos os casos
            pivo = (-Qcarter - (lz-a*E)**2 - m*m*r*r + (util*util/delta) )/delta
            print("pivo",pivo)
            if((pivo<0) and (np.abs(pivo)>10**(-30))):
                print("ERRO: Trajetória proibida R(r)<0 ",pivo)
                print("Condição Inicial será p_r = 0")
                self.x[4] = 0.0 
            else:
                self.x[4] = np.sign(self.x[4]) * np.sqrt(pivo) #pr
             
            pivo = Qcarter + (lz-a*E)**2 - (((lz-a*E*seno*seno)/seno)**2) - m*m*a*a*cosseno*cosseno
            if((pivo<0) and (np.abs(pivo)>10**(-30)) ):
                print("ERRO: Trajetória proibida Theta(theat)<0 ",pivo)
                print("Condição inicial será ptheta = 0")
                self.x[5] = 0.0
            else:
                self.x[5] = np.sign(self.x[5]) * np.sqrt(np.abs(pivo)) #ptheta 
        
        # Adiciona a dinâmica que rege essa linha de mundo 
        if(dudt == None):
            # Caso a partícula não receba uma dinâmica específica, utilizar a dinâmica padrão de órbitas em Kerr-Newman
            self.dudt = Linha_de_Mundo.dudt0 
        else:
            self.dudt = dudt
        
        
             
    def evolve_RK(self): # Método de Runge-Kutta 4 para frente no tempo
        self.xt.append(self.x) # Adiciona estado dinâmico à lista da história da partícula 
        x,y,z = self.coordinates()
        self.xplot.append(x)
        self.yplot.append(y)
        self.zplot.append(z)
        #x -- Variáveis dinâmicas (coordendas e momentos conjugados) x = ( t, r , theta , phi , p_r , p^theta )
        
        #OBS: As variáveis dinâmicas t, e phi não são exatamente variáveis dinâmicas, no entanto, é necessário atualizà-las
        
        #const -- Valor das constantes de movimento [ Energia, Momento angular em torno de z,  Massa , Carga elétrica ] 
        #dt -- Passo temporal da discretização temporal do método,permite plasticidade
        
        k1 = self.dudt(self.x,self.const) 

        k2 = self.dudt((self.x)+(k1*(self.dt)/2),self.const) 
        
        k3 = self.dudt((self.x)+(k2*(self.dt)/2),self.const)
        
        k4 = self.dudt((self.x)+(k3*(self.dt)),self.const)
        
        
        un = (k1+2*k2+2*k3+k4)/6 # Quadri-velocidade média estimada 
        
        
        self.x = (self.x) + (self.dt)*un
        
    
    
    def evolve_RK_reverse(self): # Método de Runge-Kutta para trás no tempo
        self.xt.append(self.x) # Adiciona estado dinâmico à lista da história da partícula 
        x,y,z = self.coordinates()
        self.xplot.append(x)
        self.yplot.append(y)
        self.zplot.append(z)
        #x -- Variáveis dinâmicas (coordendas e momentos conjugados) x = ( t, r , theta , phi , pr , utheta )
        
        #OBS: As variáveis dinâmicas t, e phi não são exatamente variáveis dinâmicas, no entanto, é necessário atualizà-las
        
        #const -- Valor das constantes de movimento [ Energia, Momento angular em torno de z,  Massa , Carga elétrica ] 
        #dt -- Passo temporal da discretização temporal do método,permite plasticidade
        
        k1 = -self.dudt(self.x,self.const) 

        k2 = -self.dudt((self.x)+(k1*(self.dt)/2),self.const) 
        
        k3 = -self.dudt((self.x)+(k2*(self.dt)/2),self.const)
        
        k4 = -self.dudt((self.x)+(k3*(self.dt)),self.const)
        
        un = (k1+2*k2+2*k3+k4)/6 # Quadri-velocidade média estimada
        
        self.x = (self.x) + (self.dt)*un
    
    
class Space_Time: # Essa classe representa o conjunto de linhas de mundo presente no universo e que evolui e plota gráficos de linhas de mundo
        
    def __init__(self):
        """Cria um espaço-tempo vazio, essa classe irá representar o conjunto de linhas-de-mudno presente no 
        universo e conseguirá evoluir todas as linhas-de-mundo que contêm e plotar gráficos e trajetórias"""
        self.__Wlines: List[Linha_de_Mundo] = []
        self.__active_Wlines: List[Linha_de_Mundo] = []
        self.__plot_Wlines: List[Linha_de_Mundo] = []
    def append(self,ln:Linha_de_Mundo):
        """ Adiciona Linha de mundo ao espaço-tempo
        Entrada: ln Linha_de_Mundo"""
        self.__Wlines.append(ln)
        self.__active_Wlines.append(ln)
        
    def deactivate_ln(self, ln:Linha_de_Mundo ):
        self.__active_Wlines.remove(ln)
    
    def remove_ln(self, ln:Linha_de_Mundo):
        """ Remove Linha_de_mundo do espaço-tempo
        Entrada: ln Linha_de_Mundo"""
        self.__Wlines.remove(ln)
        self.__active_Wlines.remove(ln)
    
    def generate_space_time(self): # Desenvolve todas as trajetórias 
        """Evolui todas as linhas-de-mundo pertecentes ao espaço-tempo, isso pode demorar um pouco"""
        for ln in self.__Wlines:
            print("Evoluindo trajetória:",self.__Wlines.index(ln)+1)
            print("rbk:",ln.r_horizon)
            if(ln.reverse == 0):
                for i in range(ln.N):
                    ln.evolve_RK()
                    
                    #Encontrou Buraco negro ou passou próximo demais para erros numéricos
                    if ln.x[1] <= ln.r_horizon * 1.001:
                       print("Partícula Capturada ou Muito próximo do horizonte de eventos") 
                       ln.N = i
                       break
                   
            else:
                for i in range(ln.N):
                    ln.evolve_RK_reverse()
    
    def plot_WLines(self,title = 'Órbitas de partículas testes em torno do Buraco Negro de Kerr-Newman',M = 0, a = 0, Q = 0,names = None , inside = None):
        """Plota as trajetórias das partículas-testes
        Entrada:
        title : String -> Título do plot das trajetórias
        names : String[] -> Lista contendo o nome de cada uma das trajetórias
        inside: Int -> Qualquer valor diferente de None colocará as legendas dentro da imagem, 
        1: Para upper_left, 2: upper right, 3: lower_left, 4: lower_right
        M: Float -> Energia do buraco negro  
        a: Float -> Momento angular por massa do buraco negro
        Q: Float -> Carga elétrica do buraco negro
        """
        
        plt.style.use('dark_background')
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim([-20,20])
        ax.set_ylim([-20,20])
        ax.set_zlim([-20,20])
        
        
        
        # Plota Linhas de mundo
        N = len(self.__Wlines)
        colors = plt.cm.gist_rainbow(np.linspace(0, 1, N, endpoint=False))
        i = 0
        for ln in self.__Wlines:
            ax.plot(ln.xplot,ln.yplot,ln.zplot, lw = 0.7, color = colors[i])
            i = i + 1
    
        theta = np.linspace(0,  np.pi, 300)
        phi = phi = np.linspace(0, 2*np.pi, 300)
        theta,phi = np.meshgrid(theta,phi)
        # Plota buraco negro
        r_bh = M + np.sqrt(M**2 - a**2 - Q**2) 
        x = r_bh * np.sin(theta) * np.cos(phi)
        y = r_bh * np.sin(theta) * np.sin(phi)
        z = r_bh * np.cos(theta)
        ax.plot_surface(x, y, z, color = "black")    
        
        
        #Plota ergosfera
        def r_erg(theta): return (M + np.sqrt(M**2 - (a*np.cos(theta))**2 - Q**2))
        
        r_ergosfera = r_erg(theta)
        x_erg = r_ergosfera * np.sin(theta) * np.cos(phi)
        y_erg = r_ergosfera * np.sin(theta) * np.sin(phi)
        z_erg = r_ergosfera * np.cos(theta) 
        ax.plot_surface(x_erg,y_erg,z_erg, color = "pink", alpha = 0.3) 
        ax.view_init(elev=90, azim=0)
        # Legenda
        if(names != None): 
            names.append("Buraco Negro")
            names.append("Ergosfera")
            if(inside == None):
                plt.legend(names, loc='upper left', bbox_to_anchor=(1.05, 1))
            else:
                if(inside == 1):
                    plt.legend(names, loc='upper left') 
                elif(inside == 2):
                    plt.legend(names,loc = 'upper right')
                elif(inside == 3):
                    plt.legend(names,loc = 'lower left')
                elif(inside == 4):
                    plt.legend(names,loc = 'lower right')        
        plt.title(title)
        plt.show()
        
    def plot_Graph(self,index):
        """ Plota gráficos de Wlines com index especificados
        Entrada:
        index: int[] -> Lista contendo índices das linhas-de-mundo em _Wlines que se quer plotar a dinâmica
        """
        
        
        fig, axs = plt.subplots(3, 1)
        N = len(self.__Wlines)
        colors = plt.cm.gist_rainbow(np.linspace(0, 1, N, endpoint=False))
        axs[0].set_title("r(t)")
        axs[1].set_title("θ(t)")
        axs[2].set_title("φ(t)")
        
        for k in index:
            lm = self.__Wlines[k]
            t = []
            r = []
            theta = []
            phi = []
            for j in range(lm.N):
                t.append(lm.xt[j][0])
                r.append(lm.xt[j][1])
                theta.append(lm.xt[j][2])
                phi.append(lm.xt[j][3])
            
            #Plotar gráficos em função de t
            axs[0].plot(t, r, 'o', color='red')  # 'o' = círculo
        
            axs[1].plot(t,theta, lw = 0.7, color = colors[k])
        
            axs[2].plot(t,phi, lw = 0.7, color = colors[k])

        plt.tight_layout()
        plt.show()
    
    def returnWline(self,i):
        return self.__Wlines[i]

        
