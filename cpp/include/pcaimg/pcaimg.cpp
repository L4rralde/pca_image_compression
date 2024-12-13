#include <math.h>
#include "pcaimg.h"
#include <stdexcept>
#include <iostream>

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
