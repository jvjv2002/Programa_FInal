import numpy as np
import matplotlib.pyplot as plt

# Carrega o arquivo
data = np.loadtxt("shadow.dat")

# Termos 0 devem ser pretos e 1 brancos
data_inv = 1 - data
# Plota como imagem

plt.imshow(data_inv, cmap="binary", origin="lower")
plt.xlabel("0,2 M")
plt.ylabel("0,2 M")
plt.colorbar()
plt.show()