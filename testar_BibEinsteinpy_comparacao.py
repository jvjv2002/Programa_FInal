
import numpy as np
import sys
print(sys.executable)
from einsteinpy.metric import KerrNewman
from einsteinpy.coordinates import BoyerLindquistDifferential
from einsteinpy.bodies import Body
from einsteinpy.geodesic import Geodesic
# parâmetros do buraco negro
M = 1.0           # definir sua unidade (ex: GM/c² = 1)
a = 0.5           # spin
Q = 0.3           # carga

# condição inicial para partícula/teste
r0 = 5
theta0 = np.pi/2 
phi0 = 0.0 
Lz = 0.009165100000000002
E = 
#[0, 5, 1.5707963267948966, 0, 0.0, 3.3424135056647803]
# velocidades iniciais ou momento — depende se você faz timelike ou null
v_r0 = 0.0      # ajuste conforme sua normalização
v_theta0 = 3.3424135056647803/(r0**2)
v_phi0 = 0.0

# definir o corpo test‐partícula
BLdiff = BoyerLindquistDifferential(r=r0, theta=theta0, phi=phi0,
                                     v_r=v_r0, v_theta=v_theta0, v_phi=v_phi0,
                                     a=a)
Attractor = Body(name="BH", mass=M, a=a, Q=Q)
Particle = Body(differential=BLdiff, parent=Attractor)

# calcular a geodésica
geod = Geodesic(body=Particle, metric=KerrNewman, end_lambda=1000.0, step_size=0.5)
traj = geod.trajectory

# agora você pode analisar "traj" — coordenadas, momento, etc.