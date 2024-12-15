/*
Author: Emmanuel A. Larralde Ortiz
Description:
    Performs img compression based on img columns PCA
    while dumping all matrix.
*/
#include "../include/pgm/pgm.h"
#include "../include/pcaimg/pcaimg.h"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char **argv){
    //Args safeguard
    if(argc != 2)
        return 0;

    //Arg parsing
    string fpath(argv[1]);
    //Loads src img.
    Img img = read_pgm(fpath);

    //computes img shape
    int rows = img.size();
    int cols = img[0].size();

    //Get means and std deviations vector from img
    //And sotred at normalizer object.
    Normalizer normalizer(img);
    //Print means vector
    cout << "Means: " << 1 << endl;
    for(int i = 0; i < cols; ++i)
        cout << normalizer.means[i] << " ";
    cout << endl;
    //Print std deviations vector.
    cout << "Stds: " << 1 << endl;
    for(int i = 0; i < cols; ++i)
        cout << normalizer.stds[i] << " ";
    cout << endl;

    //Normalize img columns using means and std deviations vector.
    vector< vector<double> > normalized = normalizer.normalize(img);
    cout << "Normalized image: " << rows << endl;
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j)
            cout << normalized[i][j] << " ";
        cout << endl;
    }

    //Computes cov matrix of normalized columns.
    vector< vector<double> > cov_mat = cov(normalized);
    cout << "Cov Matrix: " << cols << endl;
    for(int i = 0; i < cols; ++i){
        for(int j = 0; j < cols; ++j)
            cout << cov_mat[i][j] << " ";
        cout << endl;
    }

    //Get top k dominant eigen pairs of the covariance matrix.
    int k = 10;
    Eigen eigen(cov_mat, k);
    //Print eigenvalues
    cout << "Eigenvalues: " << 1 << endl;
    for(int eigvec = 0; eigvec < k; ++eigvec)
        cout << eigen.eigen_values[eigvec] << " ";
    cout << endl;
    //Print eigenvectors
    cout << "Eigenvectors: " << k << endl;
    for(int i = 0; i < k; ++i){
        for(int j = 0; j < cols; ++j)
            cout << eigen.eigen_vectors[i][j] << " ";
        cout << endl;
    }

    //Projects normalized columns to orthormal base of eigenvectors.
    vector< vector<double> > projected = eigen.project(normalized);
    cout << "Projected: " << rows << endl;
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < k; ++j)
            cout << projected[i][j] << " ";
        cout << endl;
    }

    //Reconstructs matrix columns from projected using base of eigenvectors.
    vector< vector<double> > reconstructed = eigen.reconstruct(projected);
    cout << "Reconstructed: " << rows << endl;
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j)
            cout << reconstructed[i][j] << " ";
        cout << endl;
    }

    //Denormalizes reconstructed using mean and std_deviations vectors.
    Img denormalized = normalizer.denormalize(reconstructed);
    cout << "Denormalized: " << rows << endl;
    for(int i = 0; i < rows; ++i){
        for(int j = 0; j < cols; ++j)
            cout << static_cast<int>(denormalized[i][j]) << " ";
        cout << endl;
    }

    return 0;
}
