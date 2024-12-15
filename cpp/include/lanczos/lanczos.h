#ifndef LANCZOS_H
#define LANCZOS_H

#include <vector>

double sinc(double x);

double lanczos_window(double x, int a);

std::vector< std::vector<unsigned char> > lanczos_resize(
    std::vector< std::vector<unsigned char> >& img,
    int nrows,
    int ncols
);

#endif