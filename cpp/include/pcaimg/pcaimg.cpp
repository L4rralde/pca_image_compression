/*
Author: Emmanuel A. Larralde Ortiz
Description:
    Functions to compress images using Principal Component Analysis.
*/

#include <math.h>
#include "pcaimg.h"
#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include "../matrices/matrices.h"

//Overall mean of Img
double mean(Pgm_Img& img){
    double acc = 0;

    int rows = img.size();
    int cols = img[0].size();
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            acc += img[i][j];
    return acc/(rows * cols);
}

//img columns/rows means
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

//Matrix columns/rows means
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

//Overall std deviation of img.
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

//img cols/rows std deviations.
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

//Null constructor of Normalizer. NOP.
Normalizer::Normalizer() {}

//Constructor from img of Normalizer.
//takes the img columns means and std deviations.
Normalizer::Normalizer(Pgm_Img& img){
    means = mean(img, 0);
    stds = std_deviation(img, 0);
}

//Normalizes the columns of img using class' mean and std deviation vectors.
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

//Ensures denormalized pixels don't overflow/underflow
unsigned char clip(double x){
    if(x < 0)
        return 0;
    if(x > 255)
        return 255;
    return static_cast<unsigned char>(x);
}

//Denormalizes a normalized matrix using mean and std deviation vectors.
Pgm_Img Normalizer::denormalize(std::vector< std::vector<double> >& normalized){
    int rows = normalized.size();
    int cols = normalized[0].size();
    Pgm_Img denormalized(
        rows, std::vector<unsigned char> (cols, 0)
    );
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            denormalized[i][j] = clip(normalized[i][j]*stds[j] + means[j]);
    return denormalized;
}

//Computes the covariance matrix of the cols of matrix.
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

//Null constructor.
Eigen::Eigen() {};

//Constructor. Computes all eigen pairs from matrix.
Eigen::Eigen(std::vector< std::vector<double> >& matrix){
    _matrix = &matrix;
    if(_matrix->size() != _matrix->at(0).size())
        throw std::invalid_argument("Matrix miust be squared");
    int n = _matrix->size();
    //Future, compute all eigenpairs.
    eigen_values = std::vector<double> (n, 0.0);
    eigen_vectors = std::vector< std::vector<double> > (
        n, std::vector<double> (n, 0.0)
    );
}

//Constructor. Computes top k dominant eigenpairs of matrix.
Eigen::Eigen(
    std::vector< std::vector<double> >& matrix,
    int k
){
    _matrix = &matrix;
    if(_matrix->size() != _matrix->at(0).size())
        throw std::invalid_argument("Matrix miust be squared");
    power(k);
}

//Power method interface.
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

//Projects mat columns to the orthonormal base of eigenvectors.
std::vector< std::vector<double> > Eigen::project(
    std::vector< std::vector<double> > &mat
){
    int rows = mat.size();
    int cols = mat[0].size();
    int k = eigen_vectors.size();
    int n = eigen_vectors[0].size();
    if(n != cols)
        throw std::invalid_argument(
            "Number of elements per eigenvector != num image cols"
        );
    
    std::vector< std::vector<double> > coordinates(
        rows, std::vector<double> (k, 0.0)
    );

    for(int row = 0; row < rows; ++row)
        for(int eigvec = 0; eigvec < k; ++eigvec)
            for(int col = 0; col < cols; ++col)
                coordinates[row][eigvec] +=
                    mat[row][col] * eigen_vectors[eigvec][col];
    return coordinates; 
}

//Reconstructs a matrix from a projection using the orthonormal base of eigenvectors.
std::vector< std::vector<double> > Eigen::reconstruct(
    std::vector< std::vector<double> > &projected
){
    int k = eigen_vectors.size();
    int n = eigen_vectors[0].size();
    int rows = projected.size();
    int projected_cols = projected[0].size();
    if(k != projected_cols)
        throw std::invalid_argument(
            "Number of eigenvectors != num of cols of projection"
        );
    
    std::vector< std::vector<double> > reconstructed(
        rows, std::vector<double> (n, 0.0)
    );

    for(int row = 0; row < rows; ++row)
        for(int col = 0; col < n; ++col)
            for(int eigvec = 0; eigvec < k; ++eigvec)
                reconstructed[row][col] +=
                    projected[row][eigvec] * eigen_vectors[eigvec][col];
    
    return reconstructed;    
}
