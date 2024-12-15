/*
Author: Emmanuel A. Larralde Ortiz
Description:
    Resizes an image using lanczos resampling interpolation.
*/
#include <iostream>
#include "../include/pgm/pgm.h"
#include "../include/lanczos/lanczos.h"

using namespace std;

int main(int argc, char **argv){
    //Agrs safeguard
    if(argc < 3)
        return 0;

    //Args parsing
    string src_path(argv[1]); //Source image path
    string dst_path(argv[2]); //Destination image path

    //Load soruce image
    Img img = read_pgm(src_path);

    //get source image size
    int rows = img.size();
    int cols = img[0].size();

    //lanczos resize.
    Img new_image = lanczos_resize(img, rows*M_PI, cols*M_PI);

    //Write back new image.
    write_pgm(dst_path, new_image);

    return 0;
}