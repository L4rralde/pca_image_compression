/*
Author: Emmanuel A. Larralde Ortiz
Description:
    Functions to compress images using Principal Component Analysis.
*/

#ifndef PCAIMG_H
#define PCAIMG_H

#include <vector>

typedef std::vector< std::vector<unsigned char> > Pgm_Img;

//Overall mean of img
double mean(Pgm_Img& img);

//Img columns/rows means.
std::vector<double> mean(Pgm_Img& img, int axis);

//Matrix columns/rows means.
std::vector<double> mean(std::vector< std::vector<double> >& img, int axis);

//Overall std deviation of img.
double std_deviation(Pgm_Img& img);

//img columns/rows std deviations.
std::vector<double> std_deviation(Pgm_Img& img, int axis);

//Class which normalizes and denormalizes an image 
//given a pair of means and std deviations vectors.
class Normalizer{
private:
public:
    std::vector<double> means; //Made public for debugging
    std::vector<double> stds; //Made public for debugging
    Normalizer();
    Normalizer(Pgm_Img& img);
    //Normalizes image columns.
    std::vector< std::vector<double> > normalize(Pgm_Img& img);
    //Denormalizes image columns.
    Pgm_Img denormalize(std::vector< std::vector<double> >& normalized);
};

//Gets the covariance matrix of the columns of matrix.
std::vector< std::vector<double> > cov(
    std::vector< std::vector<double> >& matrix
);

//Class which computes eigenvectors and eigenvalues,
//projects all columns of a matrix to the orthonormal base of eigenvectors
//and reconstructs matrix columns from the projection.
class Eigen{
private:
    //Pointer to the matrix to compute eigen pairs.
    std::vector< std::vector<double> >* _matrix;
    void power(int k); //Power method
public:
    std::vector<double> eigen_values;
    std::vector< std::vector<double> > eigen_vectors;
    //Null constructor
    Eigen();
    //FUTURE. Computes all eigenvectors and eigenvalues from matrix.
    Eigen(std::vector< std::vector<double> >& matrix);
    //Top k dominant eigen pairs from matrix.
    Eigen(
        std::vector< std::vector<double> >& matrix,
        int k
    );
    //Projects mat columnst to orthonormal base of eigenvectors.
    std::vector< std::vector<double> > project(
        std::vector< std::vector<double> > &mat
    );
    //Reconstructs columns from the projection.
    std::vector< std::vector<double> > reconstruct(
        std::vector< std::vector<double> > &projected
    );
};

#endif