/*
Author: Emmanuel A. Larralde Ortiz
Description:
    Functions for lanczso image resampling.
*/

#ifndef LANCZOS_H
#define LANCZOS_H

#include <vector>

//sinc = sin(x)/x
double sinc(double x);

//L_a(x). Interpolating function for lanczos resampling.
double lanczos_window(double x, int a);

//Lanczos resampling.
std::vector< std::vector<unsigned char> > lanczos_resize(
    std::vector< std::vector<unsigned char> >& img,
    int nrows,
    int ncols
);

#endif