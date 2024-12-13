#ifndef PCAIMG_H
#define PCAIMG_H

#include <vector>

typedef std::vector< std::vector<unsigned char> > Pgm_Img;

float mean(Pgm_Img& img);
std::vector<float> mean(Pgm_Img& img, int axis);
float std_deviation(Pgm_Img& img);
std::vector<float> std_deviation(Pgm_Img& img, int axis);

class Normalizer{
private:
    std::vector<float> means;
    std::vector<float> stds;
public:
    Normalizer();
    Normalizer(Pgm_Img& img);
    std::vector< std::vector<float> > normalize(Pgm_Img& img);
    Pgm_Img denormalize(std::vector< std::vector<float> >& normalized);
};


#endif