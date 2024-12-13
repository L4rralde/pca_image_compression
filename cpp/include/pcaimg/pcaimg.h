#ifndef PCAIMG_H
#define PCAIMG_H

#include <vector>

typedef std::vector< std::vector<unsigned char> > Pgm_Img;

float mean(Pgm_Img& img);
std::vector<float> mean(Pgm_Img& img, int axis);
float std_deviation(Pgm_Img& img);
std::vector<float> std_deviation(Pgm_Img& img, int axis);


#endif