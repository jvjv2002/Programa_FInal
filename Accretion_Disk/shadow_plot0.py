import numpy as np
import matplotlib.pyplot as plt

# Carrega o arquivo
data = np.loadtxt("disco5.dat")

# Termos 0 devem ser pretos e 1 brancos
data_inv = 1 - data
# Plota como imagem

plt.imshow(data_inv, cmap="binary", origin="lower")

plt.show()