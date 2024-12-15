/*
Author: Emmanuel A. Larralde Ortiz
Description:
    Compresses an image using Principal Component Analysis
    and enhances its resolution via lanczos resampling.
*/

#include "../include/pgm/pgm.h"
#include "../include/pcaimg/pcaimg.h"
#include "../include/lanczos/lanczos.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char **argv){
    //Args safeguard
    if(argc < 3)
        return 0;
    //Agrs parsing
    string src_path(argv[1]); //Source image path
    string dst_path(argv[2]); //Destination image path
    int k = argc >3?  atoi(argv[3]): 64; //number of top dominant eigen pairs.

    //Load soruce image.
    Img img = read_pgm(src_path);

    //source image size
    int rows = img.size();
    int cols = img[0].size();

    //Calculates img columns' means and standard deviations.
    Normalizer normalizer(img);

    //Normalize img columns.
    vector< vector<double> > normalized = normalizer.normalize(img);

    //Computes covariance matrix of normalized columns.
    vector< vector<double> > cov_mat = cov(normalized);

    //Computes top k eigen pairs of the covariance matrix.
    Eigen eigen(cov_mat, k);

    //Projects normalized image to base out of eigenvectors.
    vector< vector<double> > projected = eigen.project(normalized);

    //Reconstructs a matrix using projection.
    vector< vector<double> > reconstructed = eigen.reconstruct(projected);

    //Denormalized reconstructed matrix using img cols menas and stds
    Img denormalized = normalizer.denormalize(reconstructed);

    //lanczos remsampling to upsacale the image.
    Img upscaled = lanczos_resize(denormalized, 4*rows, 4*cols);

    //Writes back the PCA-compressed-and-reconstructed image
    string fout_path = "test.pgm";
    write_pgm(fout_path, denormalized);

    //Writes back the enhanced compressed image.
    write_pgm(dst_path, upscaled);

    return 0;
}
