#ifndef PCAIMG_H
#define PCAIMG_H

#include <vector>

typedef std::vector< std::vector<unsigned char> > Pgm_Img;

double mean(Pgm_Img& img);
std::vector<double> mean(Pgm_Img& img, int axis);
std::vector<double> mean(std::vector< std::vector<double> >& img, int axis);
double std_deviation(Pgm_Img& img);
std::vector<double> std_deviation(Pgm_Img& img, int axis);

class Normalizer{
private:
    std::vector<double> means;
    std::vector<double> stds;
public:
    Normalizer();
    Normalizer(Pgm_Img& img);
    std::vector< std::vector<double> > normalize(Pgm_Img& img);
    Pgm_Img denormalize(std::vector< std::vector<double> >& normalized);
};

std::vector< std::vector<double> > cov(
    std::vector< std::vector<double> >& matrix
);

class Eigen{
private:
    std::vector< std::vector<double> >* _matrix;
public:
    std::vector<double> eigen_values;
    std::vector< std::vector<double> > eigen_vectors;
    Eigen();
    Eigen(std::vector< std::vector<double> >& matrix);

    void jacobi(); //Symmetric matrices only
    void power(int k);
};

#endif