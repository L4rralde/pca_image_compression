//matrices.c
/*
Autor: Emmanuel Alejandro Larralde Ortiz    | emmanuel.larralde@cimat.mx
Descripcion:
    Libreria con funciones varias para manejar matrices y vectores, 
    y resolver sistemas de ecuaciones.
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "matrices.h"


//Lee un vector almacenado en un archivo txt.
double * vec_from_txt(char *fname, int *len){
    //fname: ruta del archivo, len: puntero para almacenar longitud del vector.
    FILE *file = fopen(fname, "r");
    if (file == NULL)
        return NULL;
    fscanf(file, "%d", len);

    double * ptr = (double *) malloc(sizeof(double) * (*len));
    for(int i=0; i<(*len); ++i)
        fscanf(file, "%lf", ptr + i);

    fclose(file);
    return ptr;
}

//Escribe una matriz en un archivo txt
int vec_to_txt(double *v, int n, char *fname){
    //*v: vector.
    //n: Longitud del vector.
    //fname: path del archivo.
    FILE *file = fopen(fname, "w");
    if (file == NULL)
        return -1; //No se pudo abrir el archivo.
    fprintf(file, "%d\n", n);
    for(int i=0; i<n-1; ++i)
        fprintf(file, "%lf\n", v[i]);
    fprintf(file, "%lf", v[n-1]);

    fclose(file);
    return 0;
}

//Lee una matriz almacenada en un archivo de texto
double * mat_from_txt(char *fname, int *m, int *n){
    //fname: ruta del archivo.
    //m: puntero al número de renglones.
    //n: puntero al número de columnas.
    FILE *file = fopen(fname, "r");
    if (file == NULL)
        return NULL;

    fscanf(file, "%d %d", m, n);

    int len = (*m) * (*n);
    double * ptr = (double *) malloc(sizeof(double) * len);
    for(int i=0; i<len; ++i)
        fscanf(file, "%lf", ptr + i);

    fclose(file);
    return ptr;
}

//Escribe una matriz en un archivo txt
int mat_to_txt(double *mat, int m, int n, char *fname){
    //*mat: matrix.
    //n: cantidad de filas de la matriz.
    //m: Cantidad de columnas.
    //fname: path del archivo.
    FILE *file = fopen(fname, "w");
    if (file == NULL)
        return -1; //No se pudo abrir el archivo.
    
    fprintf(file, "%d ", m);
    fprintf(file, "%d\n", n);
    for(int i=0; i<m-1; ++i){
        for(int j=0; j<n-1; ++j){
            fprintf(file, "%lf ", *(mat + i*n + j));
        }
        fprintf(file, "%lf\n", *(mat + i*n + n-1));
    }
    for(int j=0; j<n; ++j)
        fprintf(file, "%lf ", *(mat + (m-1)*n + j));

    fclose(file);
    return 0;
}

//Imprime el contenido de una matriz.
void print_matrix(double *mat, int m, int n){
    //mat: puntero a la matriz, m: número de renglones, n: número de columnas.
    for(int i=0; i<m; ++i){
        for(int j=0; j<n; ++j)
            printf("%.10lf ", *(mat + i*n + j));
        printf("\n");
    }
}

//Genera una nueva matriz (memoria dinámica) con el mismo contenido que una matriz dada mat
double *copy_matrix(double *mat, int m, int n){
    int len = m * n;
    double *copy;
    copy = (double *) malloc(len * sizeof(double));
    for(int i=0; i<len; ++i)
        *(copy +i) = *(mat +i);
    return copy;
}

//Copia los elementos de una matriz mat a una matrix copy.
void copy_matrix_no_alloc(double *mat, double *copy, int m, int n){
    int len = m*n;
    for(int i=0; i<len; ++i)
        *(copy + i) = *(mat + i);
}

//Multiplica dos matrices, retorna un puntero a la matriz resultante.
double *matmul(double * mat_a, double * mat_b, int m, int n, int o){
    //mat_a: puntero a la matriz de la izquierda (en la multiplicación).
    //mat_b: puntero a la matriz de la derecha (en la multiplicación).
    //m: número de renglones de la primera matriz.
    //n: número de columnas de la primera y de renglones de la segunda.
    //o: número de columnas de la segunda matriz.
    int max = m * o;
    double *mat_res = (double *) malloc(sizeof(double) * max);
    
    for(int i=0; i < max; ++i)
        *(mat_res + i) = 0.0;
    
    for(int i=0; i<m; ++i)
        for(int j=0; j<o; ++j)
            for(int k=0; k<n; ++k)
                *(mat_res + i*o + j) += (*(mat_a + i*n + k))*(*(mat_b + k*o + j));
    
    return mat_res;
}

//Multiplica dos matrices y guarda el resultado en otra.
void matmul_noalloc(double *mat_res, double * mat_a, double * mat_b, int m, int n, int o){
    //mat_res: puntero a la matriz resultado.
    //mat_a: puntero a la matriz de la izquierda (en la multiplicación).
    //mat_b: puntero a la matriz de la derecha (en la multiplicación).
    //m: número de renglones de la primera matriz.
    //n: número de columnas de la primera y de renglones de la segunda.
    //o: número de columnas de la segunda matriz.
    int no = n * o;
    for(int i=0; i<no; ++i)
        *(mat_res + i) = 0.0;
    for(int i=0; i<m; ++i)
        for(int j=0; j<o; ++j)
            for(int k=0; k<n; ++k)
                *(mat_res + i*o + j) += (*(mat_a + i*n + k))*(*(mat_b + k*o + j));    
}

//Resuelve un sistema de ecuaciones de la forma Dx = b. Con D diagonal.
double * solve_d(double *mat, double *b, int n){
    //mat: puntero a la matriz diagonal D.
    //b: puntero al vector b.
    //n: número de renglones y columnas de la matriz, y longitud del vector.
    double *x = (double *) malloc(sizeof(double) *  n);
    for(int i=0; i<n; ++i)
        *(x + i) = (*(b + i))/(*(mat + i*(n + 1)));
    return x;
}


//Resuelve un sistema de ecuaciones de la forma Ux = b. Con U triangular superior.
double * solve_u(double *matrix, double *vector, int n){
    //martrix: puntero a la matriz triangular superior.
    //vector: puntero al vector b.
    //n: número de renglones y columnas de la matriz, y longitud del vector.
    double *x;
    x = (double *) malloc(sizeof(double) * n);
    for(int i=0; i<n; ++i)
        *(x + i) = 0.0;

    double sum;
    for(int i=n-1; i>=0; --i){
        sum = 0.0;
        for(int j=i+1; j<n; ++j)
            sum += *(matrix + i*n + j) * *(x + j);
        *(x + i) = (*(vector + i) - sum)/(*(matrix + i*(n + 1)));
    }
    return x;
}


//Resuelve un sistema de ecuaciones de la forma Lx = b. Con L triangular inferior.
double * solve_l(double *matrix, double *vector, int n){
    //martrix: puntero a la matriz triangular inferior.
    //vector: puntero al vector b.
    //n: número de renglones y columnas de la matriz, y longitud del vector.
    double *x;
    x = (double *) malloc(sizeof(double) * n);
    for(int i=0; i<n; ++i)
        *(x + i) = 0.0;

    double sum;
    for(int i=0; i<n; ++i){
        sum = 0.0;
        for(int j=0; j<i; ++j)
            sum += *(matrix + i*n +j) * *(x + j);
        *(x + i) = (*(vector + i) - sum)/(*(matrix + i*(n + 1)));
    }
    return x;
}

//Resuelve Ax = b. Con A cuadrada usando el método de eliminación Gaussiana sin pivoteo.
double * gauss(double *matrix, double *vector, int n){
    //martrix: puntero a la matriz cuadrada.
    //vector: puntero al vector b.
    //n: número de renglones y columnas de la matriz, y longitud del vector.
    double *x;
    double *mcopy = copy_matrix(matrix, n, n);
    double *vcopy = copy_matrix(vector, 1, n);

    double factor;
    for(int col=0; col<n-1; ++col){
        for(int row = col+1; row<n; row++){
            factor = (*(mcopy + row*n + col))/(*(mcopy + col*n + col));
            for(int j=col; j<n; j++)
                *(mcopy + row*n + j) -= (*(mcopy + col*n + j)) * factor;
            *(vcopy + row) -= factor * (*(vcopy + col));
        }
    }

    x = solve_u(mcopy, vcopy, n);

    free(vcopy);
    free(mcopy);
    return x;
}

//Encuentra los índices del elemento más grande de la submatriz cuadrada M[start: n][start: n].
void _find_max_pivot(double *m, int n, int start, int *i_max, int *j_max){
    //m: puntero a la matriz m.
    //n: número de renglones y columnas de la matriz, y longitud del vector.
    //start: indice del renglon y columna que delimita la submatriz M[start: n][start: n].
    //i_max: puntero al renglón donde está el elemento de magnitud más grande en la submatriz.
    //j_max: puntero a la columna donde está el elemento de magnitud más grande en la submatriz.
    int li_max, lj_max;
    li_max = start;
    lj_max = start;
    double max, current;
    max = fabs(*(m + start*(n + 1))); //Inicializa máximo en m[start][start]
    for(int i=start; i<n; ++i){
        for(int j=start; j<n; ++j){
            current = fabs(*(m + i*n + j));
            if(current > max){ //Si m[i][j] > max anterior. Actualiza max.
                max = current;
                li_max = i;
                lj_max = j;
            }
        }
    }
    //Guarda el resultado.
    *i_max = li_max;
    *j_max = lj_max;
}

//Realiza intercambios de renglones y columnas al sistema de ecuaciones Ax = b
//para que el elemento a[src_i][src_j] intercambie posición con a[dst_i][dst_j]
void swap_matrix_elements(
    double *m, int *xi, double *b,
    int n, int src_i, int src_j, int dst_i, int dst_j
){
    //m: puntero a la matriz.
    //xi: puntero al arreglo que registra el órden del vector original x.
    //b: puntero al vector b.
    //n: número de renglones y columnas de la matriz, y longitud del vector.
    //src_i: renglón de uno de uno de los elementos a intercambiar.
    //src_j: columna de uno de uno de los elementos a intercambiar.
    //dst_i: renglón del otro elemento a intercambiar.
    //dst_j: columna del otro elemento a intercambiar.
    double aux;
    //Intercambio de renglones
    for(int j=0; j<n; ++j){
        aux = *(m + dst_i*n +j);
        *(m + dst_i*n +j) = *(m + src_i*n +j);
        *(m + src_i*n +j) = aux;
    }
    aux = *(b + dst_i);
    *(b + dst_i) = *(b + src_i);
    *(b + src_i) = aux;
    
    //Intercambio de columnas
    for(int i=0; i<n; ++i){
        aux = *(m + i*n + dst_j);
        *(m + i*n + dst_j) = *(m + i*n + src_j);
        *(m + i*n + src_j) = aux;
    }
    //xi registra todos los movimientos que se efectúan en el vector x.
    int auxi;
    auxi = *(xi + dst_j);
    *(xi + dst_j) = *(xi + src_j);
    *(xi + src_j) = auxi;
}

//Método de Gauss con pivote completo
double * gauss_pivot(double *matrix, double *vector, int n){
    //martrix: puntero a la matriz cuadrada.
    //vector: puntero al vector b.
    //n: número de renglones y columnas de la matriz, y longitud del vector.
    double *x;
    double *mcopy = copy_matrix(matrix, n, n);
    double *vcopy = copy_matrix(vector, 1, n);

    int i_max, j_max;
    //xids registra todos los movimientos que se efectúan en el vector x.
    int xids[n];
    for(int i=0; i<n; ++i)
        *(xids + i) = i;

    double factor;
    for(int col=0; col<n-1; ++col){
        //Encuentra el elemento más grande de la submatrix inferior al renglón col.
        _find_max_pivot(mcopy, n, col, &i_max, &j_max);
        //Intercambia el valor más grande con el pivote actual.
        swap_matrix_elements(mcopy, xids, vcopy, n, col, col, i_max, j_max);
        //Gauss simple
        for(int row = col+1; row<n; row++){
            factor = (*(mcopy + row*n + col))/(*(mcopy + col*n + col));
            for(int j=col; j<n; j++)
                *(mcopy + row*n + j) -= (*(mcopy + col*n + j)) * factor;
            *(vcopy + row) -= factor * (*(vcopy + col));
        }
    }
    //Resulve la matriz triangular superior que quedó como resultado de la eliminación.
    x = solve_u(mcopy, vcopy, n);

    //Reorganiza x.
    double *xcopy;
    xcopy = copy_matrix(x, 1, n);
    for(int i=0; i<n; ++i)
        *(x + xids[i]) = *(xcopy + i);
    free(xcopy);

    free(vcopy);
    free(mcopy);
    return x;
}

//Calcula el cuadrado de la distancia entre dos vectores.
double distance(double *vec_a, double *vec_b, int n){
    //vec_a: puntero de un vector.
    //vec_b: puntero del otro vector.
    //n: longitud de los punteros
    double sum = 0.0;
    double diff;
    for(int i=0; i<n; ++i){
        diff = *(vec_a + i) - *(vec_b + i);
        sum += diff * diff;
    }

    return sqrt(sum);
}

//Calcula la norma euclideana de un vector.
double norm(double *v, int n){
    //v: vector
    //n: longitud del vector
    double acc = 0;
    for(int i=0; i<n; ++i)
        acc += (*(v + i)) * (*(v + i));
    return sqrt(acc);
}

//Factoriza una matriz cuadrada A nxn en dos matrices L y U por el método de Crout.
int lu_crout(double *A, double *L, double *U, int n){
    //A: matriz nxn a fatorizar.
    //L: Matriz triangular inferior resultado de la factorización.
    //U: Matriz triangular superior resultado de la factorización.
    //n: dimensión de la matriz cuadrada.
    double acc;
    *U = 1; //U[0][0] = 1 
    *L = *A; //L[0][0] = A[0][0]
    for(int i=1; i<n; ++i){
        //Resuelve L[0:i-1][0:i-1]U[0:i-1][i] = A[0:i-1][i]
        for(int ui=0; ui<i; ++ui){
            if(*(L + ui*n + ui) == 0)
                return -1;
            acc = 0;
            for(int k=0; k<ui; ++k)
                acc += (*(L + ui*n + k)) * (*(U + k*n +i));
            *(U + ui*n + i) = (*(A + ui*n + i) - acc)/(*(L + ui*n + ui));
        }
        //U[i][i] = 1;
        *(U + i*n + i) = 1.0;
        //Resuelve U[1:i][1:i]^TL[i][1:i]^T = A[i][1:i]^T
        for(int j=0; j<=i; ++j){
            acc = 0;
            for(int k=0; k<j; ++k)
                acc += (*(L + i*n + k)) * (*(U + k*n +j));
            *(L +i*n + j) = *(A +i*n + j) - acc;
        }
    }
    return 0;
}

//Factoriza una matriz simétrica cuadrada nxn definida positiva en LL' = A.
int cholesky(double *mat, double *l, int n){
    //mat: matriz a factorizar.
    //l: matriz triangular inerior talque ll' = mat
    //n: dimensión de la matriz cuadrada
    double acc, aux;
    for(int j=0; j<n; ++j){ //Por cada columna
        //Cálculo de L[j][j]
        acc = 0;
        for(int k=0; k<j; k++)
            acc += (*(l + j*n +k)) * (*(l + j*n +k));
        aux = *(mat + j*n + j) - acc;
        if(aux <= 0) //Nan safeguard.
            return -1;
        *(l + j*n + j) = sqrt(aux); //L[j][j]
        for(int i=j+1; i<n; ++i){ //Para cada elemento de la columna debajo de la diagonal.
            //Cálculo de L[i][j], i>j
            acc = 0;
            for(int k=0; k<j; ++k)
                acc += (*(l + j*n + k)) * (*(l + i*n + k));
            *(l + i*n + j) = (*(mat + j*n + i) - acc)/(*(l + j*n + j));
        }
    }
    return 0;
}

//Factoriza una matriz simétrica cuadrada nxn definida positiva mat en mat = ldl'
int cholesky_ldl(double *mat, double *l, double *d, int n){
    //mat: matriz a factorizar.
    //l: matriz triangular inerior resultado de la factorización.
    //d: matriz diagonal resultado de la factorización.
    //n: dimensión de la matriz.
    double acc, aux;
    for(int j=0; j<n; ++j){ //Por cada columna
        //Cálculo de D[j][j]
        acc = 0;
        for(int k=0; k<j; ++k)
            acc += (*(d + k*n +k)) * (*(l + j*n +k))*(*(l + j*n +k));
        aux = *(mat + j*n + j) - acc;
        if(aux == 0) //Nan safeguard
            return -1;
        *(d + j*n + j) = aux; //D[j][j]
        *(l + j*n + j) = 1; //L[j][j]
        for(int i=j+1; i<n; ++i){ //Para cada elemento de la columna y debajo de la diagonal.
            //Calculo de L[i][j], i>j
            acc = 0;
            for(int k=0; k<j; ++k)
                acc += (*(d + k*n + k)) * (*(l + i*n + k)) * (*(l + j*n + k));
            *(l +i*n + j) = (*(mat + i*n +j) - acc)/aux;
        }
    }
    return 0;
}

//Calcula la transpuesta de una matriz mxn
double *transpose(double *mat, int m, int n){
    //mat: matriz.
    //m: numero de renglones
    //n: numero de columnas
    double *mat2;
    mat2 = (double *) malloc(m * n * sizeof(double));
    for(int i=0; i<m; ++i)
        for(int j=0; j<n; ++j)
            *(mat2 + j*m + i) = *(mat + i*n +j); //mat2[j][i] = mat[i][j]
    return mat2;
}

//Realiza el metodo de Gauss-Seidel usando los valores de xs mas recientes 
//o de la estimacion anterior
double *_gauss_seidel(double *mat, double *b, int n, int no_stale){
    double *x;
    x = (double *) calloc(n, sizeof(double));
    //Safeguard: non-0s check
    for(int i=0; i<n; ++i){
        if(*(mat + i*n + i) == 0){
            printf("EXCEPTION: found zeros in diagonal\n");
            return x;
        }
    }
    double *prev_x, *xptr, acc, error;
    prev_x = (double *) calloc(n, sizeof(double));
    //if no_stale==1:
    //  Usa los valores más recientes de x
    //else:
    //  Usa los valores de la iteración anterior.
    xptr = (no_stale == 1)? x: prev_x;
    for(int rep=0; rep<200000; rep++){
        copy_matrix_no_alloc(x, prev_x, 1, n);
        for(int i=0; i<n; ++i){//Para cada elemento de x.
            acc = 0;
            for(int j=0; j<n; ++j){//Para cada elemento de x...
                if(i==j) //Excepto xi (o sea, el mismo elemento de x)
                    continue;
                acc += (*(mat + i*n + j)) * (*(xptr + j)); //mat[i][j]*mat[j]
            }
            //x[i] = (b[i] - acc)/(mat[i][i])
            *(x + i) = (*(b + i) - acc)/(*(mat + i*n +i));
        }
        //Calcula la distancia entre la iteración actual y la anterior
        error = distance(x, prev_x, n); 
        if(error <= 0.000001) //Si no cambia mucho, termina.
            break;
    }
    if(error > 0.0001)
        printf("WARNING: Couldn't find a solution\n");

    free(prev_x);
    return x;
}

//Obtiene la solucion de Ax=b usando el método de Gauss-Seidel
double *gauss_seidel(double *mat, double *b, int n){
    return _gauss_seidel(mat, b, n, 1);
}

//Obtiene la solucion de Ax=b usando el método de Jacobi
double *jacobi(double *mat, double *b, int n){
    return _gauss_seidel(mat, b, n, 0);
}

//Genera una matriz nxm con todos los elementos igual a v
double *full(int n, int m, double v){
    int len = n * m;
    double *mat = (double *) malloc(len * sizeof(double));
    for(int i=0; i<len; ++i)
        *(mat + i) = v;

    return mat;
}

//Crea una matriz idéntidad nxn
double *eye(int n){
    //n: dimensión de la matriz
    double *mat = (double *) malloc(n * n * sizeof(double));
    for(int i=0; i<n; ++i)
        for(int j=0; j<n; ++j)
            *(mat + i*n + j) = (i == j);
    return mat;
}

//Crea una matrix nxm de zeros
double *zeros(int n, int m){
    return (double *)calloc(n * m, sizeof(double));
}

//Libera la memoria solicitada para almacenar una matriz
void free_matrix(double *mat){
    free(mat);
}

//Normaliza los elementos de un vector
void normalize_vector(double *vector, int n){
    double acc = 0;
    for(int i=0; i<n; ++i)
        acc += vector[i] * vector[i];

    if(acc == 0)
        return;

    double norm = sqrt(acc);
    for(int i=0; i<n; ++i)
        vector[i] /= norm;
}

//Resuleve sistemas LUx = b
double *solve_lu(double *l, double *u, double *b, int n){
    double *ux = solve_l(l, b, n);
    double *x = solve_u(u, ux, n);
    free(ux);
    return x;
}

//Realiza el producto punto entre dos vectores
double dot(double *v1, double *v2, int n){
    double acc = 0;
    for(int i=0; i<n; ++i)
        acc += v1[i] * v2[i];
    return acc;
}

//Elimina todas las componentes en un vector de una lista de de vectores.
//Al final normaliza el vector.
void ortho_normalize(double *v, int n, double **vs, int k){
    //V: vector a modificar.
    //n: longitud de los vectores.
    //vs: lista de vectores cuyas proyecciones serán eliminadas de v.
    //k: Cantidad de vectores en vs.
    if(k == 0){
        normalize_vector(v, n);
        return;
    }
    double c;
    for(int i=0; i<k; ++i){
        c = dot(v, vs[i], n);
        for(int j=0; j<n; ++j)
            v[j] -= c*vs[i][j];
        normalize_vector(v, n);
    }
}

//Normaliza todas las columnas de una matriz m x n
void normalize_columns(double *mat, int m, int n){
    //mat: apuntador de la matriz
    //m: cantidad de filas de la matriz.
    //n: cantidad de columnas.
    double norm;
    for(int c=0; c<n; ++c){
        norm = 0.0;
        for(int r=0; r<m;  ++r)
            norm += *((mat + r*n + c)) * (*(mat + r*n + c));
        norm = sqrt(norm);
        for(int r=0; r<m;  ++r)
            *(mat + r*n + c) /= norm;
    }
}

//Encuentra los índices del elemento con magnitud más grande de una matriz.
void find_dominant(double *mat, int m, int n, int *i, int *j){
    //mat: matriz
    //m: número de filas de mat
    //n: número de columnas de mat.
    //i, j: punteros con la ubicación del elemento dominante.
    double curr;
    double max = -INFINITY;
    int max_i, max_j;
    max_i = -1;
    max_j = -1;
    for(int i=0; i<m; ++i){
        for(int j=0; j<n; ++j){
            curr = fabs(*(mat + i*n + j));
            if(curr <= max)
                continue;
            max_i = i;
            max_j = j;
            max = curr;
        }
    }
    *i = max_i;
    *j = max_j;
}

//Encuentra los índices del elemento con magnitud más grande y que no está en la diagonal.
void find_dominant_no_diagonal(double *mat, int m, int n, int *i, int *j){
    //mat: matriz
    //m: número de filas de mat
    //n: número de columnas de mat.
    //i, j: punteros con la ubicación del elemento dominante.
    double curr;
    double max = -INFINITY;
    int max_i, max_j;
    max_i = -1;
    max_j = -1;
    for(int i=0; i<m; ++i){
        for(int j=0; j<n; ++j){
            if(i == j)
                continue;
            curr = fabs(*(mat + i*n + j));
            if(curr <= max)
                continue;
            max_i = i;
            max_j = j;
            max = curr;
        }
    }
    *i = max_i;
    *j = max_j;
}

//Aplica una rotación por la izquierda en un sólo eje
void rotate_left(double *matrix, int m, int n, double th, int p, int q){
    //matrix: Matriz a rotar.
    //m: número de filas de matrix.
    //n: número de columnas de matrix.
    //th: angulo en radianes de la rotación.
    //p, q: indican el eje para hacer la rotación.
    double *copy = copy_matrix(matrix, m, n);
    double c, s;
    c = cos(th);
    s = sin(th);
    //Multiply [0, ..., c, ...,  -s0, ..., 0] with all collumns
    for(int j=0; j<n; ++j)
        *(matrix + p*n + j) = c * (*(copy + p*n + j)) - s * (*(copy + q*n + j));
    //Multiply [0, ..., s0, ...,  c0, ..., 0] with all collumns
    for(int j=0; j<n; ++j)
        *(matrix + q*n + j) = s * (*(copy + p*n + j)) + c * (*(copy + q*n +j));
    free(copy);
};

//Aplica una rotación por la derecha en un sólo eje
void rotate_right(double *matrix, int m, int n, double th, int p, int q){
    //matrix: Matriz a rotar.
    //m: número de filas de matrix.
    //n: número de columnas de matrix.
    //th: angulo en radianes de la rotación.
    //p, q: indican el eje para hacer la rotación.
    double *copy = copy_matrix(matrix, m, n);
    double c, s;
    c = cos(th);
    s = sin(th);

    for(int i=0; i<m; ++i)
        *(matrix + i*n + p) = c * (*(copy +i*n + p)) + s * (*(copy +i*n + q));
    for(int i=0; i<m; ++i)
        *(matrix + i*n + q) = -s * (*(copy +i*n + p)) + c * (*(copy +i*n + q));

    free(copy);
};

//Método de Jacobi para una factorización de una matriz cuadrada usando sus valores propios.
void jacobi_eigen(double *mat, double *e_vecs, double *e_values, double *e_vecs_t, int n){
    //mat: matriz cuadrada a factorizar
    //e_vecs: matriz de vectores propios.
    //e_values: matriz diagonal de valores propios.
    //e_veca_t: Matriz transpuesta de vectores propios.
    //n: dimensiones de la matriz cuadrada.
    const double tol = 0.0000001;
    copy_matrix_no_alloc(mat, e_values, n, n);
    double dominant, th;
    int i, j;
    //Se inicializan e_vecs y e_vecs_t como matrices identidades.
    for(int i=0; i<n; ++i){
        for(int j=0; j<n; ++j){
            *(e_vecs + i*n + j) = i == j;
            *(e_vecs_t + i*n + j) = i == j;
        }
    }
    dominant = INFINITY;
    //Repetir hasta que el elemento dominante fuera de la diagonal sea practicamente 0.
    while(fabs(dominant) > tol){
        find_dominant_no_diagonal(e_values, n, n, &i, &j);
        dominant = (*(e_values +i*n +j));
        th = 0.5 * atan2(
            2 * dominant,
            *(e_values + i*(n+1)) - *(e_values + j*(n+1))
        );
        rotate_right(e_values, n, n, th, i, j);
        rotate_right(e_vecs, n, n, th, i, j);
        rotate_left(e_values, n, n, -th, i, j);
        rotate_left(e_vecs_t, n, n, -th, i, j);
    }
}


//Determina si una matriz cuadrada es idéntidad o no.
int is_identity(double *mat, int n){
    //mat: matriz a evaluar.
    //n: dimensiones de la matriz.
    double *I = eye(n);
    double dist = distance(mat, I, n*n);
    free(I);
    return (dist < 0.0001);
}

//FIXME: Un poco feo?
//Método iterativo de la potencia.
double **potencia(double *mat, int n, int k, int reps, double tol, double *ls, double *init){
    double *v0, *v1, l1, l0, error;

    //Vector con los eigenvectores encontrados.
    double **vs = (double **) malloc(k * sizeof(double *));
    //printf("Eigenvalores:\n");
    for(int i=0; i<k; ++i){//Repetir hasta que se encuentren k eigenvalores.
        //Inicializacion para la estimacion
        error = INFINITY;
        l0 = INFINITY;
        //v0 = [1/sqrt(n), 1/sqrt(n), ..., 1/sqrt(n)]
        if(init == NULL){
            v0 = full(1, n, 1.0/sqrt(n));
        }else{
            v0 = vector_column(init, n, k, i);
        }
        for(int rep=0; rep<reps; ++rep){
            v1 = matmul(mat, v0, n, n, 1); //v1 = Av0
            l1 = dot(v1, v1, n)/dot(v0, v1, n); //Estimacion eigenvalor
            error = fabs(l1 - l0); //Diferencia entre estimaciones
            l0 = l1;
            free(v0);
            v0 = copy_matrix(v1, 1, n); //v0 <- v1
            free(v1);
            //v0 se hace ortogonal a todas los vectores propios encontrados.
            //Y normalizado
            ortho_normalize(v0, n, vs, i);
            if(error <= tol)//Si la estimacion converge, termina.
                break;
        }
        //if(error > tol){
        //    printf("No se pudo encontrar un eigenvalor\n");
        //    break;
        //}
        //Se añade eigenvector recién calculado a la lista de eigenvectores.
        vs[i] = copy_matrix(v0, n, 1);
        ls[i] = l0;
        //printf("%lf ", l0); //Imprime eigenvalor recién encontrado.
        free(v0);
    }
    //printf("\n");

    return vs;
}

//Calcula la matriz inversa de una matriz cuadrada n x n no singular.
double *inverse(double *mat, int n){
    //mat: matrix cuadrada
    //n: dimensiones de la matrizs
    int n2 = n * n;
    double *l, *u;
    l = (double *) calloc(n2, sizeof(double));
    u = (double *) calloc(n2, sizeof(double));

    //Descomponemos en LU la matrix.
    int check = lu_crout(mat, l, u, n);
    if(check == -1){
        printf("No se puede realizar la factorizacion LU\n");
        return (double *) NULL;
    }

    //Resolvemos un sistema de ecuaciones n veces.
    double *inverse = (double *) calloc(n2, sizeof(double));
    double *b = (double *) calloc(n, sizeof(double));
    double *x;
    for(int i = 0; i < n; i++){
        *(b+i) = 1.0;            

        x = solve_lu(l, u, b, n);
        for(int j=0; j<n; ++j)
            *(inverse + i*n + j) = *(x + j);
        free(x);

        *(b+i) = 0.0;
    }

    free(l);
    free(u);
    free(b);
    return inverse; 
}

//FIXME: Un poco feo?
//Método de la potencia inversa
double **potencia_inversa(double *mat, int n, int k, int reps, double tol, double *ls, double *init){
    double *v0, *v1, l1, l0, error;

    double *inv = inverse(mat, n);
    //Vector con los eigenvectores encontrados.
    double **vs = (double **) malloc(k * sizeof(double *));
    printf("Eigenvalores:\n");
    for(int i=0; i<k; ++i){//Repetir hasta que se encuentren k eigenvalores.
        //Inicializacion para la estimacion
        error = INFINITY;
        l0 = INFINITY;
        //v0 = [1/sqrt(n), 1/sqrt(n), ..., 1/sqrt(n)]
        if(init == NULL || i == 0){
            v0 = full(1, n, 1.0/sqrt(n));
        }else{
            v0 = vector_column(init, n, k, i - 1);
        }
        for(int rep=0; rep<reps; ++rep){
            v1 = matmul(inv, v0, n, n, 1); //v1 = Av0
            l1 = dot(v0, v1, n)/dot(v1, v1, n); //Estimacion eigenvalor
            error = fabs(l1 - l0); //Diferencia entre estimaciones
            l0 = l1;
            free(v0);
            v0 = copy_matrix(v1, 1, n); //v0 <- v1
            free(v1);
            //v0 se hace ortogonal a todas los vectores propios encontrados.
            //Y normalizado
            ortho_normalize(v0, n, vs, i);
            if(error <= tol)//Si la estimacion converge, termina.
                break;
        }
        //if(error > tol){
        //    printf("No se pudo encontrar un eigenvalor\n");
        //    break;
        //}
        //Se añade eigenvector recién calculado a la lista de eigenvectores.
        vs[i] = copy_matrix(v0, n, 1);
        ls[i] = l0;
        printf("%lf ", l0); //Imprime eigenvalor recién encontrado.
        free(v0);
    }
    printf("\n");
    free(inv);

    return vs;
}

//Crea un vector a partir de una columna de una matriz.
double *vector_column(double *mat, int m, int n, int c){
    //mat: Apuntador a la Matriz
    //m: numero de filas de la matriz.
    //n: numero de columnas
    //c: columna a extraer.
    double *v = (double *) calloc(m, sizeof(double));
    for(int r=0; r<m; ++r)
        v[r] = *(mat + r*n + c);
    return v;
}

//Multiplicacion por un escalar de una matriz sin alojar nueva memoria
void scale_noalloc(double s, double *mat, int m, int n){
    //s: escalar.
    //mat: matriz
    //m, n: cantidad de filas y columnas de la matriz
    int len = m * n;
    for(int i=0; i<len; ++i)
        *(mat + i) *= s;
}

//Suma acc y mat. Guarda el resultado en mat. Si subs == 1, entonces acc <- acc - mat
void sum_acc(double *acc, double *mat, int m, int n, int subs){
    int len = m * n;
    if(subs == 1){
        for(int i=0; i<len; ++i)
            *(acc + i) -= *(mat + i);
    }else{
        for(int i=0; i<len; ++i)
            *(acc + i) += *(mat + i);
    }
}

//Factorización QR de una matriz cuadrada.
int qr(double *mat, int n, double *q, double *r){
    //Mat: Matriz a factorizar.
    //n: dimensiones de la matriz.
    //q: puntero a la memoria de la matriz Q.
    //r: Puntero a la matriz R.
    double *q_col, aux_r, *a_col, *a_col_prime, norm_a_prime;
    for(int col=0; col<n; ++col){//Para cada columna
        a_col = vector_column(mat, n, n, col); //a_col: columna #col de mat.
        a_col_prime = copy_matrix(a_col, n, 1); //â_c = a_c - ...
        for(int row=0; row<col; ++row){//Desde el primer renglón hasta col-1
            q_col = vector_column(q, n, n, row); //q_col: columna #col de q.
            //r[row][col] = <q_col, a_col>
            aux_r = dot(q_col, a_col, n);
            *(r + row*n + col) = aux_r;
            // â <- â - r[row][col] q_col
            scale_noalloc(aux_r, q_col, n, 1);
            sum_acc(a_col_prime, q_col, n, 1, 1);
            free(q_col);
        }
        //r[col][col] = ||â||
        norm_a_prime = norm(a_col_prime, n);
        if(norm_a_prime == 0)
            return -1;
        *(r + col*n + col) = norm_a_prime;
        //q_[:][col] = â/||â||
        scale_noalloc(1.0/norm_a_prime, a_col_prime, n, 1);
        for(int i=0; i<n; ++i)
            *(q + i*n + col) = *(a_col_prime + i);

        free(a_col);
        free(a_col_prime);
    }
    return 0;
}


//Factorización QR de una matriz cuadrada.
int qr2(double *mat, int n, double *q, double *r){
    //Mat: Matriz a factorizar.
    //n: dimensiones de la matriz.
    //q: puntero a la memoria de la matriz Q.
    //r: Puntero a la matriz R.
    double *q_col, aux_r, *a_col, *a_col_prime, norm_a_prime;
    for(int col=0; col<n; ++col){//Para cada columna
        a_col = vector_column(mat, n, n, col); //a_col: columna #col de mat.
        a_col_prime = copy_matrix(a_col, n, 1); //â_c = a_c - ...
        for(int row=0; row<col; ++row){//Desde el primer renglón hasta col-1
            q_col = vector_column(q, n, n, row); //q_col: columna #col de q.
            //r[row][col] = <q_col, a_col>
            aux_r = dot(q_col, a_col, n);
            *(r + row*n + col) = aux_r;
            // â <- â - r[row][col] q_col
            scale_noalloc(aux_r, q_col, n, 1);
            sum_acc(a_col_prime, q_col, n, 1, 1);
            free(q_col);
        }
        //r[col][col] = ||â||
        norm_a_prime = norm(a_col_prime, n);
        if(norm_a_prime == 0)
            return -1;
        //*(r + col*n + col) = norm_a_prime;
        q_col = vector_column(q, n, n, col);
        *(r + col*n + col) = dot(q_col, a_col, n);
        free(q_col);
        //q_[:][col] = â/||â||
        scale_noalloc(1.0/norm_a_prime, a_col_prime, n, 1);
        for(int i=0; i<n; ++i)
            *(q + i*n + col) = *(a_col_prime + i);

        free(a_col);
        free(a_col_prime);
    }
    return 0;
}

//Calcula las derivadas parciales de un vector de funciones con respecto a xi
double _partial(
    double (*f)(double *),
    double *x,
    unsigned int xi
){
    const double STEP = 0.01;

    double fxH[4];
    x[xi] -= STEP;
    fxH[0] = f(x);
    x[xi] += 2*STEP;
    fxH[1] = f(x);
    x[xi] -= STEP;
    return (1/(2.0*STEP)) * (fxH[1]-fxH[0]);    
}

//Cálculo el Jacobiano de una función vectorial.
void jacobian(
    double (*f[])(double *),
    double *x,
    int flen,
    double *JM
){
    for(int i = 0; i < flen; i++)
        for(int j = 0; j < flen; j++)
            JM[i*flen + j] = _partial(f[i], x, j);
}

//Evalua una función vectorial
void eval_no_alloc(double (*f[])(double *), double *x, double *fx, int xlen){
    for(int i = 0; i < xlen; i++)
        fx[i] = f[i](x);
}

//Método de Newton para resolver sistemas de ecuaciones no lineales.
double *newton_method(
    double (*f[])(double *),
    double *x0,
    unsigned int xlen,
    double TOL,
    unsigned int flen,
    unsigned int MAX_ITERS
){
    double *x = zeros(xlen, 1); //Vector que actualizaremos
    double *fx = zeros(xlen, 1); //Vector donde guardamos la evaluacion de la funcion
    double *y = zeros(xlen, 1);  //Vector donde se guarda el resultado de resolver el sistema
    double *JM = zeros(flen, flen); //Matriz donde se guardará la matriz Jacobiana
    copy_matrix_no_alloc(x0, x, xlen, 1); //x <- x0
    for(unsigned int k = 0; k < MAX_ITERS; ++k){ //Repite MAX_ITERS veces o hasta que se encuentre una solución
        //Fx <- f(x)
        eval_no_alloc(f, x, fx, xlen);

        //Jm = Jacobian(f, x)
        jacobian(f, x, xlen, JM);

        //Resolvemos el sistema JM y = -f(x)
        for(unsigned int s = 0; s < xlen; s++) // fx <- -fx
            fx[s]*=-1;

        //y: JM y = -fx
        double *temp = gauss_seidel(JM, fx, xlen);
        copy_matrix_no_alloc(temp, y, xlen, 1);
        free(temp);

        // x <- x + y
        for(unsigned int s = 0; s < xlen; s++)
            x[s] += y[s];

        //Si ||y|| < TOL: termina
        if(norm(y, xlen) < TOL)
            break;

        //y <- x
        copy_matrix_no_alloc(x, y, xlen, 1);
    }

    //Libera la memoria dinámica usada por la función, excepto el resultado.
    free(fx);
    free(y);
    free(JM);

    return x;
}
