#include <iostream>
#include "../include/pgm/pgm.h"
#include "../include/lanczos/lanczos.h"

using namespace std;

int main(int argc, char **argv){
    if(argc < 3)
        return 0;

    string src_path(argv[1]);
    string dst_path(argv[2]);

    Img img = read_pgm(src_path);
    
    int rows = img.size();
    int cols = img[0].size();

    Img new_image = lanczos_resize(img, rows*M_PI, cols*M_PI);
    write_pgm(dst_path, new_image);

    return 0;
}