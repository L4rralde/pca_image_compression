//matrices.h
/*
Autor: Emmanuel Alejandro Larralde Ortiz    | emmanuel.larralde@cimat.mx
Descripcion:
    Header de una librería con funciones varias para manejar matrices y vectores, 
    y resolver sistemas de ecuaciones.
*/
#ifndef MATRICES_H
#define MATRICES_H

double *vec_from_txt(char *, int *);
int vec_to_txt(double *, int, char *);
double *mat_from_txt(char *, int *, int *);
int mat_to_txt(double *, int, int, char *);
void print_matrix(double *, int, int);

double *matmul(double *, double *, int, int, int);
void matmul_noalloc(double *, double *, double *, int, int, int);
double *solve_d(double *, double *, int);
double *solve_u(double *, double *, int);
double *solve_l(double *, double *, int);
double * gauss(double *, double *, int);
void swap_matrix_elements(double *, int *, double *, int, int, int, int, int);
double * gauss_pivot(double *, double *, int);

double distance(double *, double *, int);
double norm(double *, int);

int lu_crout(double *, double *, double *, int);
int cholesky(double *, double *, int);
int cholesky_ldl(double *, double *, double *, int);

double *transpose(double *, int, int);

double *gauss_seidel(double *, double *, int);
double *jacobi(double *, double *, int);

double *full(int, int, double);
double *eye(int);
double *zeros(int, int);
void free_matrix(double *);

void normalize_vector(double *, int);
double *solve_lu(double *, double *, double *, int);

double dot(double *, double*, int);
double *copy_matrix(double *, int, int);
void copy_matrix_no_alloc(double *, double *, int, int);
void ortho_normalize(double *, int, double **, int);
void normalize_columns(double *, int, int);

void find_dominant(double *, int , int, int *, int *);
void find_dominant_no_diagonal(double *, int , int, int *, int *);
void rotate_left(double *, int, int, double, int, int);
void rotate_right(double *, int, int, double, int, int);
void jacobi_eigen(double *, double *, double *, double *, int);
double **potencia(double *, int, int, int, double, double *, double *);


int is_identity(double *, int);
double *inverse(double *, int);
double **potencia_inversa(double *, int, int, int, double, double *, double *);

double *vector_column(double *, int, int, int);
void scale_noalloc(double, double *, int, int);
void sum_acc(double *, double *, int, int, int);

int qr(double *, int, double *, double *);
int qr2(double *, int, double *, double *);

double *newton_method(
    double (*f[])(double *),
    double *,
    unsigned int,
    double,
    unsigned int,
    unsigned int
);

#endif