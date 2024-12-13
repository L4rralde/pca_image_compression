#include <math.h>
#include "pcaimg.h"
#include <stdexcept>

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