#include <math.h>
#include "pcaimg.h"
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include "../matrices/matrices.h"

double mean(Pgm_Img& img){
    double acc;

    int rows = img.size();
    int cols = img[0].size();
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            acc += img[i][j];
    return acc/(rows * cols);
}

std::vector<double> mean(Pgm_Img& img, int axis){
    int rows = img.size();
    int cols = img[0].size();
    if(axis == 0){
        std::vector<double> means (cols, 0.0);
        for(int j = 0; j < cols; ++j){
            for(int i = 0; i < rows; ++i)
                means[j] += img[i][j];
            means[j] /= rows;
        }
        return means;
    }
    if(axis == 1){
        std::vector<double> means (rows, 0.0);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j)
                means[i] += img[i][j];
            means[i] /= cols;
        }
        return means;
    }
    throw std::invalid_argument("Axis must be either 0 or 1");
}

std::vector<double> mean(std::vector< std::vector<double> >& img, int axis){
    int rows = img.size();
    int cols = img[0].size();
    if(axis == 0){
        std::vector<double> means (cols, 0.0);
        for(int j = 0; j < cols; ++j){
            for(int i = 0; i < rows; ++i)
                means[j] += img[i][j];
            means[j] /= rows;
        }
        return means;
    }
    if(axis == 1){
        std::vector<double> means (rows, 0.0);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j)
                means[i] += img[i][j];
            means[i] /= cols;
        }
        return means;
    }
    throw std::invalid_argument("Axis must be either 0 or 1");
}

double std_deviation(Pgm_Img& img){
    double acc = 0;
    int rows = img.size();
    int cols = img[0].size();
    int N = rows * cols;

    double img_mean = mean(img);
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            acc += (img[i][j] - img_mean) * (img[i][j] - img_mean);

    return sqrt(acc/N);            
}

std::vector<double> std_deviation(Pgm_Img& img, int axis){
    int cols = img[0].size();
    int rows = img.size();
    std::vector<double> means = mean(img, axis);
    if(axis == 0){
        std::vector<double> stds (cols, 0.0);
        for(int j = 0; j < cols; ++j){
            for(int i = 0; i < rows; ++i)
                stds[j] += (img[i][j] - means[j]) * (img[i][j] - means[j]);
            stds[j] = sqrt(stds[j]/rows);
        }
        return stds;
    }
    if(axis == 1){
        std::vector<double> stds (rows, 0.0);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j)
                stds[i] += (img[i][j] - means[i]) * (img[i][j] - means[i]);
            stds[i] = sqrt(stds[i]/cols);
        }
        return stds;
    }
    throw std::invalid_argument("Axis must be either 0 or 1");
}

Normalizer::Normalizer() {}

Normalizer::Normalizer(Pgm_Img& img){
    means = mean(img, 0);
    stds = std_deviation(img, 0);
}

std::vector< std::vector<double> > Normalizer::normalize(Pgm_Img& img){
    int rows = img.size();
    int cols = img[0].size();
    std::vector< std::vector<double> >  normalized (
        rows, std::vector<double> (cols, 0.0)
    );
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            normalized[i][j] = (img[i][j] - means[j])/stds[j];
    
    return normalized;
}

Pgm_Img Normalizer::denormalize(std::vector< std::vector<double> >& normalized){
    int rows = normalized.size();
    int cols = normalized[0].size();
    Pgm_Img denormalized(
        rows, std::vector<unsigned char> (cols, 0)
    );
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            denormalized[i][j] =  normalized[i][j]*stds[j] + means[j];
    return denormalized;
}

std::vector< std::vector<double> > cov(
    std::vector< std::vector<double> >& matrix
){
    int rows = matrix.size();
    int cols = matrix[0].size();
    std::vector<double> means = mean(matrix, 0);

    std::vector< std::vector<double> > cov_mat(
        cols, std::vector<double> (cols, 0.0)
    );

    for(int i = 0; i < cols; ++i){
        for(int j = 0; j < cols; ++j){
            for(int k = 0; k < rows; ++k)
                cov_mat[i][j] += matrix[k][i] * matrix[k][j];
            cov_mat[i][j] /= rows - 1;
        }
    }
    return cov_mat;
}

Eigen::Eigen() {};

Eigen::Eigen(std::vector< std::vector<double> >& matrix){
    _matrix = &matrix;
    if(_matrix->size() != _matrix->at(0).size())
        throw std::invalid_argument("Matrix miust be squared");
    int n = _matrix->size();
    eigen_values = std::vector<double> (n, 0.0);
    eigen_vectors = std::vector< std::vector<double> > (
        n, std::vector<double> (n, 0.0)
    );
}

void Eigen::jacobi(){
    int n = _matrix->size();
    int n2 = n * n;

    double *matrix = (double *) calloc(n2, sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            matrix[i*n + j] = _matrix->at(i).at(j);
    double *e_vals = (double *) calloc(n, sizeof(double));
    double *e_vecs = (double *) calloc(n2, sizeof(double));
    double *e_vecs_t = (double *) calloc(n2, sizeof(double));

    std::cout << "Jacobi for " << n << std::endl;
    jacobi_eigen(matrix, e_vecs, e_vals, e_vecs_t, n);
    for(int i = 0; i < n; ++i){
        eigen_values[i] = e_vals[i];
        for(int j = 0;  j < n; ++j)
            eigen_vectors[i][j] = e_vecs[i*n + j];
    }

    free(e_vals);
    free(matrix);
    free(e_vecs);
    free(e_vecs_t);
}

void Eigen::power(int k){
    int n = _matrix->size();
    int n2 = n * n;
    double tol = 0.000001;
    int reps = 10000;
    double *ls = (double *) calloc(k, sizeof(double));
    double *matrix = (double *) calloc(n2, sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            matrix[i*n + j] = _matrix->at(i).at(j);

    double** vectors = potencia(matrix, n, k, reps, tol, ls, NULL);
    eigen_values = std::vector<double> (k, 0.0);
    eigen_vectors = std::vector< std::vector<double> > (
        k, std::vector<double> (n, 0.0)
    );
    for(int i = 0; i < k; ++i)
        eigen_values[i] = ls[i];
    for(int i = 0; i < k; ++i)
        for(int j = 0; j < n; ++j)
            eigen_vectors[i][j] =  vectors[i][j];

    free(matrix);
    free(ls);
    for(int i = 0; i < k; ++i)
        free(vectors[i]);
    free(vectors);
}