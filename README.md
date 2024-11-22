# Image compression with dominant eigenvectors and splines

1. Obtener un conjunto de datos de imágenes de mediana resolución a escala de grises en formato pgm
2. PCA:
    1. Pasar la matriz a vector plano.
    2. Eliminar la media al vector.
    3. Calcular la matriz de covarianza del vector.
    4. Calcular los $n$ eigenvectores de los $n$ eigenvalores dominantes.
    5. Proyectar el vector original al espacio de los eigenvalores dominantes.
3. Splines:
    1. Recupera la forma de la imagen.
    2. Selecciona un conjunto de puntos igualmente epspaciados para hacer una malla y solo guarda esos puntos.
4. Presentación de la imagen:
    1. Aplica la interpolación con splines para mostrar la imagen.
