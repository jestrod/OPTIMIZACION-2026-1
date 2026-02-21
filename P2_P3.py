import numpy as np
import itertools as it
import pandas as pd

# Crear matrices y vectores Ax = b
A = np.array([
        [1., 2., 4., 1., 0., 0.],
        [2., 5., 6., 0., 1., 0.],
        [1., 3., 4., 0., 0., 1.],
    ], dtype=np.float64)

x = ["x1", "x2", "x3", "s1", "s2", "s3"]

b = np.array([10., 30., 50.], dtype=np.float64)

# Arreglo para guardar las soluciones
solv = []

# Tamaño de la matrix
# n = restricciones
# m = variables
n, m = A.shape

# Iterar sobre todas las combinaciones posibles
for base in it.combinations(range(m), n):
    # Index tiene una tupla de los índices de las columnas a utilizar
    A_base = A[:, base]

    # Verificar el rango de la matriz
    # (Para saber si es invertible y saber si tiene solución)
    rank = np.linalg.matrix_rank(A_base)

    # No hay solución (no es de rango completo)
    if rank < n:
        continue

    # Resolver la ecuación
    x_base = np.linalg.solve(A_base, b)

    # Obtener vector con todas las ecuaciones
    x_full = np.zeros(m)
    x_full[list(base)] = x_base

    # Verificar si es factible
    f = np.all(x_full >= 0)

    solv.append({
        "base": [x[i] for i in base],
        **{x[i] : x_full[i] for i in range(m)},
        "factible": "Sí" if f else "No",
    })

# Crear una tabla
df = pd.DataFrame(solv)
pd.set_option('display.precision', 2)
print(df)
