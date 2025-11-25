import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import BlackHolePO as bh 
import numpy as np

sp = bh.Space_Time()

sp.plot_WLines(title = "Buraco Negro de Kerr-Newman",M = 1.0, a = 0.85, Q = 0)