"""Gera dados.txt: estatura (cm) e peso (kg) sinteticos de 1000 pessoas.

Cada linha: "<estatura> <peso>" separados por espaco, 1 casa decimal.
Semente fixa (42) para que o arquivo seja sempre o mesmo.

Uso: python gerar_dados.py   (requer numpy)
"""
import numpy as np

N = 1000
rng = np.random.default_rng(42)
estatura = rng.normal(170, 9, N)                                # ~ N(170; 9) cm
peso = 70 + 0.9 * (estatura - 170) + rng.normal(0, 9, N)        # correlacionado
estatura = np.clip(estatura, 145, 198)
peso = np.clip(peso, 42, 120)

with open("dados.txt", "w", newline="\n") as f:
    for a, p in zip(estatura, peso):
        f.write(f"{a:.1f} {p:.1f}\n")
print(f"dados.txt gerado com {N} pessoas")
