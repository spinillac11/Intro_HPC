import numpy as np
import matplotlib.pyplot as plt

# Parámetros constantes
W = 1
L = 1
F = 1

# Leer la matriz desde el archivo de texto
matrix_file_path = "matrix.txt"
matrix1 = np.loadtxt(matrix_file_path)

# Leer el vector desde el archivo de texto
vector_file_path = "vec.txt"
vector = np.loadtxt(vector_file_path)

# Seleccionar las partes relevantes de la matriz y el vector
A = matrix1[1:, 1:]  # Ignorar la primera fila y columna
b = vector[1:]      # Ignorar el primer elemento del vector

# Resolver el sistema de ecuaciones lineales
x = np.linalg.solve(A, b)

# Agregar un valor de 0 al principio del vector x
x = np.insert(x, 0, 0)

# Crear los valores de x para graficar
n = len(x)
x_values = np.linspace(0, 0.5, num=n)

# Graficar los elementos del vector resultante en función de los múltiplos de (1/2)/(n-1)
plt.plot(x_values, x, marker='o', linestyle='-', label='Aprox elementos finitos')

# Calcular y graficar la función
func_values = -(W/L)*np.power(x_values,4)/24 + (W+F)*np.power(x_values,3)/12 - (W*np.power(L,2)/24 + F*np.power(L,2)/16)*x_values
plt.plot(x_values, func_values, linestyle='--', label='Solucion analitica')

plt.xlabel('$x$')
plt.ylabel('$\eta (x)$')
plt.title('Superficie neutra')
plt.legend()
plt.grid(True)
plt.show()
