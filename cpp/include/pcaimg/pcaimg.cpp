#include <math.h>
#include "pcaimg.h"
#include <stdexcept>
#include <iostream>

float mean(Pgm_Img& img){
    float acc;

    int rows = img.size();
    int cols = img[0].size();
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            acc += img[i][j];
    return acc/(rows * cols);
}

std::vector<float> mean(Pgm_Img& img, int axis){
    int rows = img.size();
    int cols = img[0].size();
    if(axis == 0){
        std::vector<float> means (cols, 0.0);
        for(int j = 0; j < cols; ++j){
            for(int i = 0; i < rows; ++i)
                means[j] += img[i][j];
            means[j] /= rows;
        }
        return means;
    }
    if(axis == 1){
        std::vector<float> means (rows, 0.0);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j)
                means[i] += img[i][j];
            means[i] /= cols;
        }
        return means;
    }
    throw std::invalid_argument("Axis must be either 0 or 1");
}

std::vector<float> mean(std::vector< std::vector<float> >& img, int axis){
    int rows = img.size();
    int cols = img[0].size();
    if(axis == 0){
        std::vector<float> means (cols, 0.0);
        for(int j = 0; j < cols; ++j){
            for(int i = 0; i < rows; ++i)
                means[j] += img[i][j];
            means[j] /= rows;
        }
        return means;
    }
    if(axis == 1){
        std::vector<float> means (rows, 0.0);
        for(int i = 0; i < rows; ++i){
            for(int j = 0; j < cols; ++j)
                means[i] += img[i][j];
            means[i] /= cols;
        }
        return means;
    }
    throw std::invalid_argument("Axis must be either 0 or 1");
}

float std_deviation(Pgm_Img& img){
    float acc = 0;
    int rows = img.size();
    int cols = img[0].size();
    int N = rows * cols;

    float img_mean = mean(img);
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            acc += (img[i][j] - img_mean) * (img[i][j] - img_mean);

    return sqrt(acc/N);            
}

std::vector<float> std_deviation(Pgm_Img& img, int axis){
    int cols = img[0].size();
    int rows = img.size();
    std::vector<float> means = mean(img, axis);
    if(axis == 0){
        std::vector<float> stds (cols, 0.0);
        for(int j = 0; j < cols; ++j){
            for(int i = 0; i < rows; ++i)
                stds[j] += (img[i][j] - means[j]) * (img[i][j] - means[j]);
            stds[j] = sqrt(stds[j]/rows);
        }
        return stds;
    }
    if(axis == 1){
        std::vector<float> stds (rows, 0.0);
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

std::vector< std::vector<float> > Normalizer::normalize(Pgm_Img& img){
    int rows = img.size();
    int cols = img[0].size();
    std::vector< std::vector<float> >  normalized (
        rows, std::vector<float> (cols, 0.0)
    );
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            normalized[i][j] = (img[i][j] - means[j])/stds[j];
    
    return normalized;
}

Pgm_Img Normalizer::denormalize(std::vector< std::vector<float> >& normalized){
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

std::vector< std::vector<float> > cov(
    std::vector< std::vector<float> >& matrix
){
    int rows = matrix.size();
    int cols = matrix[0].size();
    std::vector<float> means = mean(matrix, 0);

    std::vector< std::vector<float> > cov_mat(
        cols, std::vector<float> (cols, 0.0)
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
