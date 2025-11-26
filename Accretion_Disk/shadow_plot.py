import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
# Carrega o arquivo
data = np.loadtxt("disco5.dat")

# Definir cores para cada valor:
# 0 → preto
# 1 → vermelho
# 2 → verde
# 3 → amarelo
# 4 → azul
colors = [
    "black",   # 0
    "red",     # 1
    "green",   # 2
    "yellow",  # 3
    "blue"     # 4
]
boundaries = [0, 1 , 2, 3, 4, 5] 
cmap = ListedColormap(colors)

norm = mcolors.BoundaryNorm(boundaries, cmap.N)
# Plot da imagem
plt.imshow(data, cmap=cmap,norm = norm ,origin="lower")
plt.colorbar(label="Cores")  # opcional
plt.title("Buraco Negro de Kerr a/M = 0.99")
plt.xlabel("x")
plt.ylabel("y")
plt.show()