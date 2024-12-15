#include "lanczos.h"
#include <math.h>
#include <iostream>

double sinc(double x){
    if(x == 0)
        return 1;
    return sin(x)/x;
}

template <typename T>
T clamp(T x, T min, T max){
    if(x < min)
        return min;
    if(x > max)
        return max;
    return x;
}

double lanczos_window(double x, int a){
    if(x == 0)
        return 1;
    if(x < -a || x > a)
        return 0;
    return sinc(M_PI * x) * sinc(M_PI*x/a);
}

std::vector< std::vector<unsigned char> > lanczos_resize(
    std::vector< std::vector<unsigned char> >& img,
    int nrows,
    int ncols
){
    int rows = img.size();
    int cols = img[0].size();

    double x_step = cols/(ncols - 1.0);
    double y_step = rows/(nrows - 1.0);

    std::vector< std::vector<unsigned char> > xresampled(
        rows, std::vector<unsigned char>(ncols, 0)
    );

    for(int row = 0; row < rows; ++row){
        double x = 0;
        for(int ncol = 0; ncol < ncols; ++ncol){
            int floor_x = x;
            double dnpixel = 0;
            for(int i = floor_x - 2; i <= floor_x + 3; ++i){
                dnpixel += img[row][clamp(i, 0, cols - 1)]*
                            lanczos_window(x - i, 3);
            }
            xresampled[row][ncol] = clamp<double>(dnpixel, 0.0, 255.0);
            x += x_step;
        }
    }

    std::vector< std::vector<unsigned char> > nimg(
        nrows, std::vector<unsigned char>(ncols, 0)
    );

    for(int ncol = 0; ncol < ncols; ++ncol){
        double y = 0;
        for(int nrow = 0; nrow < nrows; ++nrow){
            int floor_y = y;
            double dnpixel = 0;
            for(int i = floor_y - 2; i <= floor_y + 3; ++i){
                dnpixel += xresampled[clamp(i, 0, rows-1)][ncol]*
                            lanczos_window(y - i, 3);
            }
            nimg[nrow][ncol] = clamp<double>(dnpixel, 0.0, 255.0);
            y+= y_step;
        }
    }

    return nimg;
}
