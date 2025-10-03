import matplotlib.pyplot as plt
import numpy as np

# Tamaño del scan (N puntos)
N = np.logspace(2, 5, 100)  # de 1e2 a 1e5 puntos

# Tiempos relativos (ficticios, solo para ilustrar)
time_inloop = N * 1.0e-3        # HBHS dentro del loop
time_outloop = N * 2.0e-4       # HBHS fuera del loop
time_post = N * 5.0e-5          # HBHS post-process

plt.figure(figsize=(8,6))

plt.loglog(N, time_inloop, label="HBHS dentro del loop", lw=2, color="red")
plt.loglog(N, time_outloop, label="HBHS fuera del loop", lw=2, color="orange")
plt.loglog(N, time_post, label="HBHS post-process", lw=2, color="green")

plt.xlabel("Número de puntos en el scan (N)")
plt.ylabel("Tiempo relativo (unidades arbitrarias)")
plt.title("Escalamiento del tiempo de evaluación con HBHS")
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("plots/hbhs_efficiency.png", dpi=150)